from roots import Root
from weyl import weyl
from basis_vectors import Vector_array
from scalar_products import scalarp_roots
from travel_diagrams import travel

class Dynkin:
    def __init__(self,v_arr,root_list):
        self.v_arr=v_arr
        self.vect_list=v_arr.vect_list
        self.root_list=root_list
        self.roots_coeffs=self.roots_coeffs()
        self.reps=self.reps()
        self.cartan_matrix=self.cartan_matrix()

    def reps(self):
        diagram=''
        roots=''
        for element in self.root_list:
            if element.length==0:
                diagram+='X'
            else:
                diagram+='    O    '
            diagram+='    '
            roots+=element.reps()+'    '
        return diagram,roots
    
    def roots_coeffs(self):
        root_coeffs=[]
        for root in self.root_list:
            root_coeffs.append(root.coeffs)
        return root_coeffs
    
    def view(self):
        print(self.reps[0])
        print(self.reps[1])

    def normalise_cartan(self,cartan_matrix):
        norm_cartan=[]
        for row in cartan_matrix:
            fractionals=[abs(el) for el in row if abs(el)<1 and el!=0]
            if fractionals:
                to_mult=1/min(fractionals)
            else:
                to_mult=1
            new_row=[int(to_mult*el) for el in row]
            norm_cartan.append(new_row)
        return norm_cartan

    def cartan_matrix(self):
        matrix=[]
        for root1 in self.root_list:
            row=[]
            for root2 in self.root_list:
                if root1.length!=0:
                    row.append(2*scalarp_roots(root1,root2)/root1.length)
                else:
                    row.append(scalarp_roots(root1,root2))
            if root1.length==0:
                row_copy=row
                min_row=min(row_copy)
                if min_row==0:
                    row_copy=[el for el in row_copy if el != 0]
                    min_row=min(row_copy)
                row=[el/min_row for el in row]
            matrix.append(row)
        return self.normalise_cartan(matrix)
    
    def print_cartan(self):
        for row in self.cartan_matrix:
            print(row)

    def weyl_on_diagram(self,root_i):
        new_root_list=[]
        if root_i>len(self.root_list):
            raise ValueError
        refl_root=self.root_list[root_i-1]
        for element in self.root_list:
            new_root_list.append(weyl(refl_root,element))
        return Dynkin(self.v_arr,new_root_list)
    
    def find_Q1(self):
        coeffs=self.roots_coeffs
        index=intersect_roots(coeffs[0],coeffs[1])
        if not index:
            index=intersect_roots(coeffs[0],coeffs[2])
        if index[0][0]==1:
            return f'Q()|({index[0][1]})'
        return f'Q ({index[0][1]})|()'

class Distinguished(Dynkin):
    def __init__(self,v_arr):
        self.v_arr=v_arr
        self.vect_list = v_arr.vect_list
        self.root_list = self.populate()
        super().__init__(self.v_arr, self.root_list)
    def populate(self):
        root_l=[]
        coeff_tuple=[1,-1]
        while len(coeff_tuple)<len(self.vect_list):
            coeff_tuple.append(0)
        for i in range(len(self.vect_list)-1):
            root_l.append(Root(tuple(coeff_tuple),self.v_arr))
            coeff_tuple.insert(0,0)
            coeff_tuple.pop()
        # last root
        coeff_tuple[-2]=1
        root_l.append(Root(tuple(coeff_tuple),self.v_arr))
        return root_l
    
#auxiliary functions
def invert(tuples):
    list=[]
    for i in range(len(tuples)):
        list.append(-tuples[i])
    return list

def intersect_roots(list1,list2):
    listinv=invert(list2)
    list=[]
    if len(list1)!=len(listinv):
        raise ValueError
    for i in range(len(list1)):
        if listinv[i]==list1[i] and listinv[i]!=0:
            list.append((i+1,list1[i]))
    return list