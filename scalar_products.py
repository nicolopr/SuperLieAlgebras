from basis_vectors import *

def scalarp_vector(vector1,vector2):
    if vector1.grading==0 & vector2.grading==0:
        if vector1.index==vector2.index:
            return 1
    elif vector1.grading==1 & vector2.grading==1:
        if vector1.index==vector2.index:
            return -1
    return 0

def scalarp_roots(root1,root2):
    res=0
    if len(root1.coeffs)==len(root2.coeffs) and root1.vect_array==root2.vect_array:
        for i in range(len(root1.coeffs)):
            res+=root1.coeffs[i]*root2.coeffs[i]*scalarp_vector(root1.vect_array.vect_list[i],root2.vect_array.vect_list[i])
        return res
    else:
        raise ValueError
    
# v_arr=Vector_array((2,1))
# root=Root((-1,3,1),v_arr)
# print('root1')
# root.view()
# root2=Root((2,1,1),v_arr)
# smroot=sum_root(root,root2)
# print('root2')

# root2.view()
# print('root1+2')

# smroot.view()
# print(scalarp_roots(root,root2))
# print(scalarp_roots(root,smroot))
# print(scalarp_roots(smroot,root2))