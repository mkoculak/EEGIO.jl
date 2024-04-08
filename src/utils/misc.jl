# Decide which file access method to use
function _read_method(fid::IO, method, type, shape)
    if method == :Direct
        return fid
    elseif method == :Mmap
        return Mmap.mmap(fid, type, shape)
    else
        error("Selected read method '$method' is not recognized. Only valid options are :Direct and :Mmap.")
    end
end

# Helper functions for decoding string and numerical header entries
_decodeString(fid, size) = strip(ascii(String(read!(fid, Array{UInt8}(undef, size)))))
_decodeNumber(fid, numType, size) = parse(numType, ascii(String(read!(fid, Array{UInt8}(undef, size)))))

# Helper function to decode channel specific string entries
function _decodeChanStrings(fid, nChannels, size)
    arr = Array{String}(undef, nChannels)
    buf = read(fid, nChannels*size)
    for i = eachindex(arr)
        arr[i] = strip(ascii(String(buf[(size*(i-1)+1):(size*i)])))
    end
    return arr
end

# Helper function to decode channel specific numerical entries
function _decodeChanNumbers(fid, numType, nChannels, size)
    arr = Array{numType}(undef, nChannels)
    buf = read(fid, nChannels*size)
    for i = eachindex(arr)
        arr[i] = parse(numType, ascii(String(buf[(size*(i-1)+1):(size*i)])))
    end
    return arr
end

# Write the general data information
function _write_record(fid, field, fieldLength; default="")
    # Prepare the entry.
    if field == ""
        record = rpad(string(default), fieldLength)
    else
        if length(field)>fieldLength
            @warn "Header field \"$field\"
                    is longer than required $fieldLength bytes and will be truncated."
        end
        record = rpad(string(field),fieldLength)
    end

    # Write bytes to file
    write(fid, codeunits(record))
end

# Write the chennel specific information
function _write_channel_records(fid, nChannels, field, fieldLength; default="")
    for chan in 1:nChannels
        #Prepare the entry for each channel.
        if field == ""
            record = rpad(string(default), fieldLength)
        else
            if length(field[chan])>fieldLength
                @warn "Header field \"$field\" entry on position $chan
                        is longer than required $fieldLength bytes and will be truncated."
            end
            record = rpad(string(field[chan]), fieldLength)
        end

        # Write bytes to file
        write(fid, codeunits(record))
    end
end

# Calculate scaling and offsets for EDF/BDF data while checking if including them makes sense.
function _resolve_offsets(header, addOffset, numPrecision)
    if addOffset & (numPrecision <: Integer)
        @warn "Reading data as integers is not compatible with scaling and offset correction,
         so it will be omitted. To get corrected data, use floating point output, e.g. NumPrecision=Float64.
         To hide this warning include addOffest=false in the read call."

        addOffset = false
    end

    if addOffset
        # Correct for misspecified physical or digital spans
        physSpan = header.physMax .- header.physMin
        replace!(x -> x <= 0 ? 1 : x, physSpan)
        digSpan = header.digMax .- header.digMin
        replace!(x -> x <= 0 ? 1 : x, physSpan)

        scaleFactors = numPrecision.(physSpan ./ digSpan)
        offsets = _get_offsets(header, numPrecision, scaleFactors)
    else
        scaleFactors = ones(numPrecision, header.nChannels)
        offsets = numPrecision.(header.physMin .* 0)
    end

    return scaleFactors, offsets
end

function _get_offsets(header, numPrecision::Type{ <: AbstractFloat}, scaleFactors)
    return numPrecision.(header.physMin .- (header.digMin .* scaleFactors))
end

function _get_offsets(header, numPrecision::Type{ <: Integer}, scaleFactors)
    return round.(numPrecision, header.physMin .- (header.digMin .* scaleFactors))
end

# Calculate buffer size close to requested memory size and number of records that fit in it.
function _get_buffer_size(memSize, recSize)
    buffSize = memSize + recSize/2
    buffSize = Int(buffSize - buffSize%recSize)
    buffSize == 0 ? buffSize = recSize : nothing
    nRecords = Int(buffSize/recSize)

    return buffSize, nRecords
end