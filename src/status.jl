import ValueHistories: push!

function checkstatus(::NoStatus, z)
    return false
end

printstatusheader(::NoStatus) = nothing
