# TEモードに関してスラブ導波路をFEMで求める
# 前提として導波路構造の媒質が等方性・非磁性(比透磁率＝1)
# 無損失（電流密度と電荷密度がJ=ρ=0）の誘電体からなるとする　　
# 光波の伝搬方向はZ,
# x方向に対して構造は不均一

# 1次元有限要素法でTEモードの電界分布を求める

module SlabWaveguide_FEM

# λ:波長
# k0:波数
# 2a:導波路幅
# n1:コア屈折率
# n2:クラッド屈折率
# n:自然数
λ = 1.0
k0 = 2pi * λ
a = 5.0
n1 = 1.49
n2 = 1.48
n = 0

# 解析解で求めたβ
#β=9.35820167262037

export All_element, make_matrix_KM

# 構造体を定義
mutable struct All_element
    N::Int64 # nodeの数
    seg::Array{Array{Float64,1},1} # 線要素を格納する配列
    K::Array{Float64,2} # 全体の剛性行列
    M::Array{Float64,2} # 全体の質量行列

    function All_element(x)
        N = length(x)
        seg = make_segments(x)
        A = zeros(Float64,N, N)
        B = zeros(Float64,N, N)
        return new(N, seg)
    end

    function make_segments(vec_x)
        N = length(vec_x)                                           # vec_xの要素数
        vec_segments = Array{Array{Float64,1},1}(undef,N-1)         # r_i,r_i+1 を配列に格納する．配列数はelementの数(vec_x-1になる)分 
        for i=1:N-1                                                 
            vec_segments[i] = [vec_x[i],vec_x[i+1]]                 
        end
        return vec_segments
    end

end

function make_matrix_KM(self)
    seg = self.seg
    N = self.N
    C = zeros(Float64,N, N)

    for e=1:N-1

        tmp_element = Element(seg[e], e, N)

        Ce = tmp_element.Ce

        C[e, e]     = Ce[1,1]
        C[e+1, e]   = Ce[2,1]
        C[e, e+1]   = Ce[1,2]
        C[e+1, e+1] = Ce[2,2]

    end

    return C

end

struct Element
    seg::Array{Float64,1}
    len::Float64
    Ce::Array{Float64,2}

    # 単一要素を記述した構造体
    function Element(seg, e, N)
        len = calc_length(seg)
        Ce = C_element(len, seg)
        
        return new(seg, len, Ce)
    end

    function calc_length(vec_segments)
        edge_n = length(vec_segments) - 1 # エッジ数=ノード数 - 1
        len = zeros(Float64,edge_n)

        # 2次元の有限要素法の場合は下記のコードになる(かも・・・)
        # for i=1:(edge_n)
        #     len[i] = vec_segments[i][2]-vec_segments[i][1]
        # end

        len = vec_segments[2]-vec_segments[1]

        return len
    end

    # 各nodeでのCeを計算する．計算対象とする構造により異なる
    function C_element(len, vec_segments)
        
        if -a ≤ vec_segments[1] &&  vec_segments[2] ≤　a  
            n_e = 1.49
        else
            n_e = 1.48
        end
        
        c11 = -6/len^2 + n_e^2 * k0
        c22 = c11
        c12 = 6/len^2
        c21 = c12        
        
        return [c11 c12; c21 c22]
    end
    
end

end

using Plots
using .SlabWaveguide_FEM
using LinearAlgebra

x_min = -10
x_max = 10
step = 0.1                                                   
x = collect(x_min:step:x_max)

ele = All_element(x)
C = make_matrix_KM(ele)
β, vec = eigen(C)

plt = plot(x, vec[:,201])
plot(plt)

#ve = eigvecs(C)
#ve = ve[201]
#β = eigmin(C)
#β = sqrt(complex(β[1]))
#println(β)




# for i in 1:length(β)
#     if β[i] > 0
#         print(i)
#         #plt = plot(r, vec[:,i])
#         #plot(plt)
#         # β_min = β[i]
#         # vec_min = vec[i]
#         break
#     end
# end 

# xx = collect(x_min:step:x_max)


# # print(β_min)
# # k0 = k0[1]
# # k0 = complex(k0)
# # k0 = sqrt(k0)
# # #println(β)
# # println(k0)

# seg = All_element.seg
# N = All_element.N
# 17.9600542967557im
# 17.963023604311534
