module yk

using Base.GMP.MPZ: sizeinbase, mpz_t
using Base.GMP: CulongMax, ClongMax, libgmp

export iroot

if Clong == Int32
    const LargeInt = Union{Int64, Int128}
    const LargeUInt = Union{UInt64, UInt128}
else
    const LargeInt = Int128
    const LargeUInt = UInt128
end



"""
    iroot(x,n)::BigInt

if `x` is a positive number,it returns the integer `n`-th root of `x`, defined as the largest integer `r` s.t. `r^n ≤ x`.

if `x` is a negative number,it returns `-iroot(-x,n)`

# Examples
```julia
julia>iroot(1331, 3)
11

julia>iroot(30, 3)
3

julia>iroot(-120, 3)
-4

julia>iroot(BigInt(2)^1000, 1000)
2

```

# Arguments

- `x`: an integer (can be negative if `n` is odd).

- `n`: a positive integer (`n > 0`).

# Notes
This function uses GMP's `__gmpz_root` for efficiency.

If `n` is extremely large and `x` is small in bit length, the result is simply `sign(x)`.
"""
function iroot(x::BigInt, n::CulongMax)::BigInt
    iszero(n) && throw(DomainError((x, n), "iroot(x, n) is undefined for n <= 0"))
    if x < 0 && iseven(n)
        throw(DomainError((x,n), "iroot(x, n) is  undefined for x<0 && even n"))
    end
    z = BigInt()
    ccall((:__gmpz_root, libgmp), Cint, (mpz_t, mpz_t, Culong), z, x, n)
    return z
end

function iroot(x::BigInt, n::ClongMax)::BigInt
    n <= 0 && throw(DomainError((x, n), "iroot(x, n) is undefined for n <= 0"))
    return iroot(x,Culong(n))
end

function iroot(x::BigInt,n::Union{LargeInt, LargeUInt})::BigInt
    if n > typemax(Culong)
        if sizeinbase(x,2) <= n
            if x < 0 && iseven(n)
                throw(DomainError((x,n), "iroot(x, n) is  undefined for x<0 && even n"))
            else
                return sign(x)
            end
        else
            throw(OverflowError("Arguments x,n passed iroot(x,n) are too large"))
        end
    end
    n <= 0 && throw(DomainError((x, n), "iroot(x, n) is undefined for n <= 0"))
    return iroot(x,Culong(n))
end

function iroot(x::BigInt,n::BigInt)::BigInt
    if iszero(ccall((:__gmpz_fits_ulong_p, libgmp), Cint, (mpz_t,), n))
        n <= 0 && throw(DomainError((x, n), "iroot(x, n) is undefined for n <= 0"))
        if sizeinbase(x,2) <= n
            if x < 0 && iseven(n)
                throw(DomainError((x,n), "iroot(x, n) is  undefined for x<0 && even n"))
            else
                return sign(x)
            end
        else
            throw(OverflowError("Arguments x,n passed iroot(x,n) are too large"))
        end
    end
    return iroot(x,Culong(n))
end

iroot(x::T, n::Integer) where {T<:Integer} = T(iroot(BigInt(x), n))

end # yk
