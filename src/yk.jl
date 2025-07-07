module yk

import Base.GMP:libgmp
import Base.GMP.MPZ.mpz_t


function iroot(x::BigInt,n::UInt32)
    x < 0 && throw(DomainError(x, "iroot(x, n) is undefined for negative x"))
    z = BigInt()
    ccall((:__gmpz_root, libgmp), Cint, (mpz_t, mpz_t, Culong), z, x, n)
    return z
end

end
