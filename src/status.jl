import ValueHistories: push!

function checkstatus(::NoStatus, z)
    return false
end

function logextra(::NoStatus, z, override=false)
    return false
end

printstatusheader(::NoStatus) = nothing
