
function write_bdf(f::String, bdf::BDF; overwrite=false, kwargs...)
    if isfile(f)
        if overwrite
            rm(f)
        else
            error("File $f already exists. Use `overwrite=true` to overwrite the file.")
        end
    end
    open(f, "w+", lock = false) do fid
        write_bdf(fid, bdf; kwargs...)
    end
end

function write_bdf(fid::IO, bdf::BDF; useOffset=true, method=:Direct, buffer=512_000, tasks=1)
    # Write the header
    newHeader = write_bdf_header(fid, bdf)

    # Write the data
    write_bdf_data(fid, bdf, newHeader, useOffset, method, buffer, tasks)
end

function write_bdf_header(fid, bdf::BDF)

    header = deepcopy(bdf.header)

    # Add Status channel data
    _add_status!(header)

    write(fid, UInt8(header.idCodeNonASCII))
    _write_record(fid, header.idCode, 7, default="BIOSEMI")
    _write_record(fid, header.subID, 80)
    _write_record(fid, header.recID, 80)
    _write_record(fid, header.startDate, 8)
    _write_record(fid, header.startTime, 8)
    _write_record(fid, header.nBytes, 8)
    _write_record(fid, header.versionDataFormat, 44, default="24BIT")
    _write_record(fid, header.nDataRecords, 8, default="-1")
    _write_record(fid, header.recordDuration, 8, default="1")
    _write_record(fid, header.nChannels, 4)
    chans = header.nChannels
    _write_channel_records(fid, chans, header.chanLabels, 16)
    _write_channel_records(fid, chans, header.transducer, 80)
    _write_channel_records(fid, chans, header.physDim, 8)
    _write_channel_records(fid, chans, header.physMin, 8)
    _write_channel_records(fid, chans, header.physMax, 8)
    _write_channel_records(fid, chans, header.digMin, 8)
    _write_channel_records(fid, chans, header.digMax, 8)
    _write_channel_records(fid, chans, header.prefilt, 80)
    _write_channel_records(fid, chans, header.nSampRec, 8)
    _write_channel_records(fid, chans, header.reserved, 32)

    return header
end

function _add_status!(header)
    header.nChannels += 1
    push!(header.chanLabels, "Status")
    push!(header.transducer, "Trigger and Status")
    push!(header.physDim, "Boolean")
    push!(header.physMin, -8388608)
    push!(header.physMax, 8388607)
    push!(header.digMin, -8388608)
    push!(header.digMax, 8388607)
    push!(header.prefilt, "No filtering")
    push!(header.nSampRec, 2048)
    push!(header.reserved, "Reserved")
end

# Write data into a file.
function write_bdf_data(fid, bdf, header, useOffset, method, buffer, tasks)

    scaleFactors, offsets = _resolve_offsets(header, useOffset, eltype(bdf.data))

    records = header.nDataRecords
    channels = header.nChannels
    srate = header.nSampRec[1]
    duration = header.recordDuration
    samples = srate * duration

    if isempty(bdf.status)
        status = Status(zeros(Int8, samples), zeros(Int8, samples), zeros(Int8, samples))
    else
        status = bdf.status
    end

    output = _read_method(fid, method, Vector{UInt8}, records * channels * samples * 3)

    write_bdf_data(output, bdf.data, status, records, channels, samples, scaleFactors, offsets, buffer, tasks)

    return nothing
end

# Write using memory mapping
function write_bdf_data(output::Vector{UInt8}, data, status, records, channels, samples, scaleFactors, offsets, buffer, tasks)
    
    @tasks for rec in 1:records
        @set ntasks = tasks

        write_record(output, data, status, rec, channels, samples, scaleFactors, offsets)
    end

    # Necessary to properly close the mmaped file.
    finalize(output)
    GC.gc(false)

    return nothing
end

function write_record(output, data, status::BDFStatus, rec, channels, samples, scaleFactors, offsets)
    rIdx = (rec - 1) * channels * samples
    dataIdx = (rec - 1) * samples

    # Skip the last channel as it is the status channel
    for chan in 1:(channels - 1)
        cIdx = (chan - 1) * samples
        for sample in 1:samples
            pointer = (rIdx + cIdx + sample - 1) * 3 + 1
            dIdx = dataIdx + sample
            recode_value!(data, output, pointer, dIdx, chan, scaleFactors, offsets)
        end
    end

    # Write the status channel
    cIdx = (channels - 1) * samples
    for sample in 1:samples
        pointer = (rIdx + cIdx + sample - 1) * 3 + 1
        dIdx = dataIdx + sample
        recode_value!(status, dIdx, output, pointer)
    end

    return nothing
end

function write_bdf_data(fid::IO, data, status, records, channels, samples, scaleFactors, offsets, buffer, tasks)
    recordSize = channels * samples * 3

    # Testing different methods of writing to disk has shown that writing in chunks,
    # is more efficient than writing the whole array at once.
    # Here we determine the size of chunk that is a multiple of record size and closest
    # to 512k which is used as a default.
    buffSize, recNum = _get_buffer_size(buffer, recordSize)
    
    output = Vector{UInt8}(undef, buffSize)

    for chunk in collect(chunks(1:records, size=recNum))
        pointer = 1
        for rec in chunk
            rIdx = (rec-1) * samples
            # Skip the last channel as it is the status channel
            for chan in 1:(channels - 1)
                for sample in 1:samples
                    dIdx = rIdx + sample
                    recode_value!(data, output, pointer, dIdx, chan, scaleFactors, offsets)
                    pointer += 3
                end
            end

            # Write the status channel
            for sample in 1:samples
                dIdx = rIdx + sample
                recode_value!(status, dIdx, output, pointer)
                pointer += 3
            end
        end
        @inbounds write(fid, view(output, 1:(length(chunk)*recordSize)))
    end

    return nothing
end

function recode_value!(data::Matrix{Int32}, output::Vector{UInt8}, pointer, dIdx, chan, scaleFactors, offsets)
    @inbounds output[pointer] = data[dIdx, chan] % UInt8
    @inbounds output[pointer+1] = (data[dIdx, chan] >> 8) % UInt8
    @inbounds output[pointer+2] = (data[dIdx, chan] >> 16) % UInt8
end

function recode_value!(data::Matrix{Int64}, output::Vector{UInt8}, pointer, dIdx, chan, scaleFactors, offsets)
    @inbounds value = round(Int32, data[dIdx, chan])
    @inbounds output[pointer] = value % UInt8
    @inbounds output[pointer+1] = (value >> 8) % UInt8
    @inbounds output[pointer+2] = (value >> 16) % UInt8
end

function recode_value!(data::Matrix{<:AbstractFloat}, output::Vector{UInt8}, pointer, dIdx, chan, scaleFactors, offsets)
    @inbounds @fastmath value = round(Int32, ((data[dIdx, chan]-offsets[chan])/scaleFactors[chan]))
    @inbounds output[pointer] = value % UInt8
    @inbounds output[pointer+1] = (value >> 8) % UInt8
    @inbounds output[pointer+2] = (value >> 16) % UInt8
end

function recode_value!(status::BDFStatus, dIdx, output::Vector{UInt8}, pointer)
    @inbounds output[pointer] = status.low[dIdx]
    @inbounds output[pointer+1] = status.high[dIdx]
    @inbounds output[pointer+2] = status.status[dIdx]
end