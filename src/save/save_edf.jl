function write_edf(f::String, edf::EDF; overwrite=false, kwargs...)
    if isfile(f)
        if overwrite
            rm(f)
        else
            error("File $f already exists. Use `overwrite=true` to overwrite the file.")
        end
    end
    open(f, "w+", lock = false) do fid
        write_edf(fid, edf; kwargs...)
    end
end

function write_edf(fid::IO, edf::EDF; useOffset=true, method=:Direct, tasks=1)
    # Write the header
    write_edf_header(fid, edf)

    write_edf_data(fid, edf.header, edf.data, useOffset, method, tasks)
end

function write_edf_header(fid, edf::EDF)
    
    header = edf.header

    _write_record(fid, header.version, 8)
    _write_record(fid, header.patientID,80)
    _write_record(fid, header.recordingID, 80)
    _write_record(fid, header.startDate, 8)
    _write_record(fid, header.startTime, 8)
    _write_record(fid, header.nBytes, 8)
    _write_record(fid, header.reserved44, 44)
    _write_record(fid, header.nDataRecords, 8)
    _write_record(fid, header.recordDuration, 8)
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
    _write_channel_records(fid, chans, header.reserved32, 32)
end

function write_edf_data(fid, header, data, useOffset, method, tasks)
    scaleFactors, offsets = _resolve_offsets(header, useOffset, Float64)

    records = header.nDataRecords
    channels = header.nChannels
    samples = header.nSampRec

    output = _read_method(fid, method, Vector{Int16}, sum(records .* samples))

    write_edf_data(output, data, records, samples, channels, scaleFactors, offsets, tasks)
end

# Write using direct IO
function write_edf_data(output::IO, data, records, samples, channels, scaleFactors, offsets, tasks)
    
    buffer = Vector{Int16}(undef, sum(samples))
    recIdx = Vector{Int}(undef, length(samples))
    recCum = 0

    for record in 1:records
        recIdx .= (record-1) .* samples
        write_record(buffer, data, recIdx, recCum, samples, channels, scaleFactors, offsets)

        write(output, buffer)
    end
end

# Multithreaded write using Mmap
function write_edf_data(output::Vector{Int16}, data, records, samples, channels, scaleFactors, offsets, tasks)
    
    @tasks for record in 1:records
        @set ntasks = tasks
        @local begin
            recIdx = Vector{Int}(undef, length(samples))
            recCum = 0
        end

        recIdx .= (record-1) .* samples
        recCum = sum(recIdx)
        write_record(output, data, recIdx, recCum, samples, channels, scaleFactors, offsets)
    end
    # Necessary to properly close the mmaped file.
    finalize(output)
end

function write_record(output, data, recIdx, recCum, samples, channels, scaleFactors, offsets)
    for chan in 1:channels
        chanIdx = sum(view(samples, 1:(chan-1)))
        for sample in 1:samples[chan]
            sampIdx = recCum + chanIdx + sample
            recode_value(data, output, recIdx, sampIdx, sample, chan, offsets, scaleFactors)
        end
    end
end

function recode_value(data::Vector{Vector{T}}, output::Vector{Int16}, recIdx, sampIdx, sample, chan, offsets, scaleFactors) where T <: AbstractFloat
    output[sampIdx] = round(Int16, (data[chan][recIdx[chan]+sample] - offsets[chan])/scaleFactors[chan])
end

function recode_value(data::Vector{Vector{T}}, output::Vector{Int16}, recIdx, sampIdx, sample, chan, offsets, scaleFactors) where T <: Integer
    output[sampIdx] = data[chan][recIdx[chan]+sample]
end