module FEM

#定数定義
p = 0 # p:TE,TMモードの場合はp=0となる
ref_ind = 1.5 # ガラスの屈折率
ε = ref_ind^2　# 比誘電率＝屈折率^2
mode = "EH" # TEモードとする(p=0ならTEモードもTMモードも同じ形になる)　 
β = 0.5 # 
l = 1
# function cal_l(p, mode)
#     if p == 0
#         l = 1
#     elseif p >= 1
#         if mode == "EH"
#             l = p +1
#         elseif mode == "HE"
#             l = p - 1
#         else
#             println("mode is undefined") 
#         end                
#     else
#         println("p is undefined")
#     end
#     return l
# end

# l = cal_l(p, mode)


    export All_element, make_matrix_KM

    # abstract type Set_constant end
    # struct constant <: Set_constant 
    #     p::Integer  
    #     ε::Float64
    #     mode::String
    #     β::Float16

    #     function constant(p, ε, mode, β)
    #         return new(p, ε, mode, β)
    #     end
    # end

    # function (constant::Type{<::Constan})(p, ε, mode, β)
    #     p = p     
    #     ε = ε
    #     mode = mode
    #     β = β
    #     return p, ε, mode, β
    # end

    # 構造体を定義
    mutable struct All_element
        N::Int64
        seg::Array{Array{Float64,1},1}
        K::Array{Float64,2}
        M::Array{Float64,2}
    
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
        K = zeros(Float64,N, N)
        M = zeros(Float64,N, N)
        for e=1:N-1
            tmp_element = Element(seg[e])

            Ke = tmp_element.Ke
            Me = tmp_element.Me

            K[e, e]     = Ke[1,1]
            K[e+1, e]   = Ke[2,1]
            K[e, e+1]   = Ke[1,2]
            K[e+1, e+1] = Ke[2,2]
                
            M[e, e]     = Me[1,1]
            M[e+1, e]   = Me[2,1]
            M[e, e+1]   = Me[1,2]
            M[e+1, e+1] = Me[2,2]
        end

        return K, M

    end

    struct Element
        seg::Array{Float64,1}
        len::Float64
        Ke::Array{Float64,2}
        Me::Array{Float64,2}

        # 単一要素を記述した構造体
        function Element(seg)
            len = calc_length(seg)
            Ke = K_element(seg, len)
            Me = M_element(seg, len)
            
            return new(seg, len, Ke, Me)
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
    
        function K_element(seg, len)
            integrals = cals_integral(seg, len)
            # All_element.lの形ではコンストラクタを呼ぶことはできない。当たり前だったけど。
            # All_elementのコンストラクタをelement構造体で呼ぶには継承する必要がある。
            return 2π * (integrals[2] + l * integrals[4] + β^2 * integrals[1])
        end
    
        function M_element(seg, len)
            integrals = cals_integral(seg, len)
            return 2π * ε * integrals[1]
        end    

        function cals_integral(segments, len)
            r1 = segments[1]
            r2 = segments[2]
             
            u = r1/len
            v = log(r2/len)
    
            A11 = -u - 1.5 + (u + 1)^2 * v
            A12 = u + 0.5 - (u + 1) * u * v
            A21 = A12
            A22 = -u + 0.5 + u^2 * v        
    
            integral1 = (r1 * len / 12) * [3 1;1 1] + (r2 * len / 12) * [1 1;1 3]
            integral2 = (r1 / 2len) * [1 -1;-1 1] + (r2 / 2len) * [1 -1;-1 1]
            integral3 = (r1 / 6) * [-2 -1;2 1] + (r2 / 6) * [-1 -2;1 2]
            integral4 = [A11 A12;A21 A22]
    
            return [integral1, integral2, integral3, integral4] 
        end

    end

end

using Plots
using .FEM
using LinearAlgebra

M = 100                                                             # 分割数
L = 90                                                            # rの最大値
r = collect(range(0,L,length=M))　                                  # r軸上の各elementの座標を格納した1次元配列

ele = All_element(r) # , p, mode, ε, β)
K, M = make_matrix_KM(ele)

k0, vec = eigen(K, M)
#k0 = sqrt(k0.values)

plt = plot(r, vec[:,1])
plot(plt)
#println(vec[:,1])


# seg = All_element.seg
# N = All_element.N

# #println(inv(felements.B)*felements.A)

# fp = open("eigvalues.dat","w")
# for i=1:length(F.values)
#     println(fp,i,"\t",F.values[i],"\t",π^2*i^2/(L)^2)
# end
# close(fp)
# fp = open("vector.dat","w")
# for i=2:length(x)-1
#     println(fp,x[i],"\t",F.vectors[i-1,1])
# end
# close(fp)
