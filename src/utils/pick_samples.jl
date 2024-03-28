"""
    EEGIO.pick_samples(header::BDFHeader, records)
    EEGIO.pick_samples(header::BDFHeader, samples, nDataSamples::Ineger)

Internal function used to properly pick the range of samples that should be read from every
channel. Since data is stored differently in each format, specialized variants were written
to read as much data as needed to satify the query without compromising efficiency.
Please refer to the description of `timeSelect` keyword argument of particular read function
for more detailed account of the picking behaviour.

Always returns a UnitRange.

## Examples

```julia
# Picking the first 10 records of the data (assuming `file` is an object of type BDF).
pick_samples(file.header, 1:10)
```
```julia
# Picking the first 10 seconds of the data (assuming `file` is an object of type EEG).
pick_samples(file.header, (1.,10.), size(file.data,1))
```
"""
function pick_samples(header, samples::Any)
    error("Selection of time interval \"$samples\" should be a number, a range, or a list of indices.")
end
 
# BDF and EDF selection
# Picking the time interval to load, measured as number of records or seconds.

function pick_samples(header::Header, records::Symbol)
    if records == :All
        return 1:_sample_count(header)
    else
        error("Unknown symbol :$records passed. Did You mean :All?")
    end
end

# Integer interpreted as an index of a data record to be read.
function pick_samples(header::Header, record::Integer)
    if 0 < record < _sample_count(header)
        return record:record
    else
        error("Number of a record to read should be between 1 and $(header.nDataRecords). Got $record instead.")
    end
end

# Unitrange interpreted as an interval including records with such indexes.
function pick_samples(header::Header, records::UnitRange)
    if records[1] >= 1 && records[end] <= _sample_count(header)
        return records
    else
        error("Range $records does not fit in the available $(_sample_count(header)) records.")
    end
end

# Tuple of floats interpreted as seconds. Picking all records which parts are included
# in the given time interval. Therefore actual data might be slightly larger then the interval.
function pick_samples(header::Header, records::Tuple{AbstractFloat, AbstractFloat})
    dur = _sample_duration(header)
    signalTime = _signal_duration(header)
    if 0. <= records[1] && records[2] <= signalTime
        return Int64(floor(records[1]/dur))+1:Int64(ceil(records[2]/dur))
    else
        error("Time range $records does not fit the available length of the data: $signalTime")
    end
end

# Range of floats should be parsed to a tuple.
pick_samples(header::Header, records::StepRangeLen) = pick_samples(header, (records[1], records[end]))

_sample_count(header::BDFHeader) = header.nDataRecords
_sample_count(header::EDFHeader) = header.nDataRecords
_sample_count(header::EEGHeader) = header.common["NumberOfSamples"]
_sample_count(header::SETHeader) = header.pnts

_sample_duration(header::BDFHeader) = header.recordDuration
_sample_duration(header::EDFHeader) = header.recordDuration
_sample_duration(header::EEGHeader) = 1_000_000 / header.common["SamplingInterval"]
_sample_duration(header::SETHeader) = 1 / header.srate

_signal_duration(header::Header) = _sample_count(header) * _sample_duration(header)
