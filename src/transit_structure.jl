mutable struct Transit_Struct{T}
  # Structure to hold arrays for transit:
  r :: T
  b :: T
  u_n :: Array{T,1}
  n :: Int64
  v_max :: Int64
  c_n :: Array{T,1}
  sn :: Array{T,1}
  Iv ::  Array{T,1}
  Jv ::  Array{T,1}
  grad:: Bool
  dIvdk ::  Array{T,1}
  dJvdk ::  Array{T,1}
  s2_grad :: Array{T,1}
  dsndr ::  Array{T,1}
  dsndb ::  Array{T,1}
  dcdu :: Array{T,2}
  dfdrbc :: Array{T,1}
  dfdrbu :: Array{T,1}
end

function transit_init(r::T,b::T,u_n::Array{T,1},grad::Bool) where {T <: Real}
# Initializs a transit structure.
n = length(u_n)
if iseven(n)
  v_max = round(Int64,n/2)+2
else
  v_max = round(Int64,(n-1)/2)+2
end
trans = Transit_Struct{T}(r,b,u_n,n,v_max,
  zeros(T,n+1),    # c_n
  zeros(T,n+1),  # sn
  zeros(T,v_max+1),# Iv
  zeros(T,v_max+1),# Jv
  true,            # grad
  zeros(T,v_max+1),# dIvdk
  zeros(T,v_max+1),# dJvdk
  zeros(T,2),      # s2_grad
  zeros(T,n+1),  # dsndr
  zeros(T,n+1),  # dsndb
  zeros(T,n+1,n),  # dcdu
  zeros(T,n+3),    # dfdrbc
  zeros(T,n+2)     # dfdrbu
)
return trans
end
