@inline function iroot(x::BigInt, n::CulongMax)::BigInt
    iszero(n) && throw(DomainError((x, n), "iroot(x, n) is undefined for n <= 0"))
    if x < 0 && iseven(n)
        throw(DomainError((x,n), "iroot(x, n) is  undefined for x<0 && even n"))
    end
    z = BigInt()
    ccall((:__gmpz_root, libgmp), Cint, (mpz_t, mpz_t, Culong), z, x, n)
    return z
end