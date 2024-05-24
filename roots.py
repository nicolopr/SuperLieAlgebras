from basis_vectors import Vector_array
from scalar_products import scalarp_roots

class Root:
    def __init__(self,coeffs: tuple, vect_array: Vector_array):
        self.coeffs=coeffs
        self.vect_array=vect_array
        self.length=scalarp_roots(self,self)
    
    def reps(self):
        res=''
        if len(self.coeffs)==len(self.vect_array.vect_list):
            for i in range(len(self.coeffs)):
                if self.coeffs[i]!=0:
                    if self.coeffs[i]>0:
                        res+='+'
                    if abs(self.coeffs[i])>1:
                        res+=f'{int(self.coeffs[i])}*{self.vect_array.vect_list[i].reps()}'
                    elif self.coeffs[i]==1:
                        res+=f'{self.vect_array.vect_list[i].reps()}'
                    elif self.coeffs[i]==-1:
                        res+=f'-{self.vect_array.vect_list[i].reps()}'
            if res:
                return res.strip('+')
            else:
                raise ValueError
        else:
            raise ValueError
    
    def view(self):
        print(self.reps())
        
def sum_root(root1: Root,root2: Root):
    if len(root1.coeffs)==len(root2.coeffs) and root1.vect_array==root2.vect_array:
        new_coeffs=[root1.coeffs[i]+root2.coeffs[i] for i in range(len(root2.coeffs))]
        return Root(new_coeffs, root1.vect_array)
    else:
        raise ValueError
    
def mult_root(root: Root, scalar: float):
    new_coeffs=[scalar*root.coeffs[i]for i in range(len(root.coeffs))]
    return Root(new_coeffs, root.vect_array)