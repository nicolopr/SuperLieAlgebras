from roots import Root
from weyl import weyl
from basis_vectors import Vector_array
from scalar_products import scalarp_roots
from travel_diagrams import travel
from copy import deepcopy

class Dynkin:
    def __init__(self,v_arr,root_list):
        self.v_arr=v_arr
        self.vect_list=v_arr.vect_list
        self.root_list=root_list
        self.roots_coeffs=self.roots_coeffs()
        self.Qindices = self.Q_indices_reorder()
        self.reps=self.reps()
        self.cartan_matrix=self.cartan_matrix()
    
    def write_indices_Q(self,list_indices,root_number):
        string='Q(' if root_number<sum(self.v_arr.bosons_fermions)-2 else 'S('
        for index in list_indices[0]:
            string+=index
        string+=')|('
        for index in list_indices[1]:
            string+=index
        return string+')'


    def roots_coeffs(self):
        root_coeffs=[]
        for root in self.root_list:
            root_coeffs.append(root.coeffs)
        return root_coeffs

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
    
    def Q_indices(self):
        unique_vectors,indices=[],[]
        bosons=self.v_arr.bosons_fermions[0]
        fermions=self.v_arr.bosons_fermions[1]
        #gl-type tail
        for i in range(bosons+fermions-3):
            unique_vectors.append(unique_elements(self.root_list[i].coeffs,self.root_list[i+1].coeffs))
        for element in unique_vectors:
            indices.append(get_indices_from_unique_vectors(bosons,fermions,indices,element))
        #fork node
        fork_node=bosons+fermions-3
        if 2 not in self.root_list[fork_node+1].coeffs and -2 not in self.root_list[fork_node+1].coeffs:
            #S+ node is not a long root
            unique_vectors.append(unique_elements(self.root_list[fork_node].coeffs,self.root_list[fork_node+1].coeffs))
        else:
            unique_vectors.append(unique_elements(self.root_list[fork_node].coeffs,self.root_list[fork_node+2].coeffs))
        indices.append(get_indices_from_unique_vectors(bosons,fermions,indices,unique_vectors[-1]))
        #spinor nodes
        if 2 not in self.root_list[fork_node+1].coeffs and -2 not in self.root_list[fork_node+1].coeffs:
            #S+ node is not a long root
            not_unique_elements_sp=not_unique_elements(self.root_list[fork_node].coeffs,self.root_list[fork_node+1].coeffs)
            unique_elements_sm=unique_elements(self.root_list[fork_node+1].coeffs,self.root_list[fork_node].coeffs)
            unique_elements_sp=[-el for el in unique_elements_sm]
            temp_indices=[get_indices_from_unique_vectors(bosons,fermions,indices,not_unique_elements_sp)]
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,unique_elements_sp))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,unique_elements_sm))
        else:
            #S+ node is a long root
            #terrible naming conventions here
            not_unique_elements_sp=not_unique_elements(self.root_list[fork_node].coeffs,self.root_list[fork_node+2].coeffs)
            unique_elements_sm=unique_elements(self.root_list[fork_node+2].coeffs,self.root_list[fork_node].coeffs)
            unique_elements_sp=[-el for el in unique_elements_sm]
            temp_indices=[get_indices_from_unique_vectors(bosons,fermions,indices,not_unique_elements_sp)]
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,unique_elements_sm))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,unique_elements_sp))
        return indices
    
    def Q_indices_reorder(self):
        Q_reordered=[]
        for Q_fun in self.Q_indices():
            new_Qfun=[]
            for set_indices in Q_fun:
                new_Qfun.append(sorted(set_indices))
            Q_reordered.append(new_Qfun)
        return Q_reordered
    
    def reps(self):
        diagram=''
        roots=''
        Q_f=''
        Qindex=self.Qindices
        for root_number,element in enumerate(self.root_list):
            if element.length==0:
                diagram+='X'
            else:
                diagram+='    O    '
            diagram+='    '
            roots+=element.reps()+'    '
            Q_f+=self.write_indices_Q(Qindex[root_number],root_number)+'  '
        return diagram,roots,Q_f
    
    def view(self):
        print(self.reps[0])
        print(self.reps[1])
        print(self.reps[2])

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
            list.append(i)
    return list

def unique_elements(list1, list2):
    list_return=list(list1).copy()
    for i in range(len(list_return)):
        if list_return[i]==-list2[i]:
            list_return[i]=0
        else:
            list_return[i]=int(list_return[i])
    return list_return

def not_unique_elements(list1, list2):
    list_return=list(list1).copy()
    for i in range(len(list_return)):
        if list_return[i]!=-list2[i]:
            list_return[i]=0
        else:
            list_return[i]=-int(list_return[i])
    return list_return

def find_nonzero(list):
    for i,element in enumerate(list):
        if element!=0:
            return i,element

def get_indices_from_unique_vectors(bosons,fermions,indices,unique_vector):
    no_vector, value_vector=find_nonzero(unique_vector)
    indices_vector=deepcopy(indices[-1]) if indices else [[],[]]
    if no_vector>=fermions:
        if value_vector>0:
            indices_vector[0].append(f'{no_vector-fermions+1}')
        elif value_vector<0:
            indices_vector[0].append(f"{no_vector-fermions+1}'")
    else:
        if value_vector>0:
            indices_vector[-1].append(f'{no_vector+1}')
        elif value_vector<0:
            indices_vector[-1].append(f"{no_vector+1}'")
    return indices_vector

# vector_array=Vector_array((3,2))
# distinguished=Distinguished(vector_array)
# visited=travel(distinguished)
# distinguished.Q_indices()
# for element in visited[:100]:
#     element.view()
