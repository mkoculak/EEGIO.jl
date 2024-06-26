module EEGIO

using Mmap
using MAT
using OhMyThreads: @tasks, @set, @local, chunks


# TODO make a test to check if this value is updated
# with Pkg.TOML.parsefile(joinpath(pkgdir(EEGIO), "Project.toml"))["version"]
version = "0.1.0"

abstract type EEGData end
abstract type Header end

# Misc functions
include("utils/misc.jl")

# BDF files
include("formats/BDF.jl")
include("load/load_bdf.jl")
include("save/save_bdf.jl")
export BDF, BDFHeader, read_bdf, write_bdf

# EDF files
include("formats/EDF.jl")
include("load/load_edf.jl")
include("save/save_edf.jl")
export EDF, EDFHeader, read_edf, write_edf

# EEG files
include("formats/EEG.jl")
include("load/load_eeg.jl")
include("save/save_eeg.jl")
export EEG, EEGHeader, EEGMarkers, read_eeg, write_eeg

# SET files
include("formats/SET.jl")
include("load/load_set.jl")
export SET, SETHeader, read_set

# Convenient type unions
BEDFHeader = Union{BDFHeader, EDFHeader}

# Helper functions
include("utils/pick_channels.jl")
include("utils/pick_samples.jl")
include("utils/resize.jl")
export crop, select, split

end # module