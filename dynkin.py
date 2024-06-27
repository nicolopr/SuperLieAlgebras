import numpy as np
from roots import Root
from weyl import weyl
from basis_vectors import Vector_array
from scalar_products import scalarp_roots
from travel_diagrams import travel
from copy import deepcopy

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
            indices_vector[0].append(f"{no_vector-fermions+1}^")
    else:
        if value_vector>0:
            indices_vector[-1].append(f'{no_vector+1}')
        elif value_vector<0:
            indices_vector[-1].append(f"{no_vector+1}^")
    return indices_vector

def reorder(root_list):
        root_list_cp=deepcopy(root_list)
        if 2 in root_list[-2].coeffs or -2 in root_list[-2].coeffs:
            #reorder long-root in sp diagrams
            temp_root=root_list_cp[-2]
            root_list_cp[-2]=root_list_cp[-1]
            root_list_cp[-1]=temp_root
        elif 2 not in root_list[-1].coeffs and -2 not in root_list[-1].coeffs:
            #need to reorder so that last node in S- node, but avoid that if diagram is sp-type
            common_coeff_to_fork=not_unique_elements(root_list[-3].coeffs,root_list[-2].coeffs)
            plus_or_minus_one=[coeff for coeff in common_coeff_to_fork if coeff!=0][0]
            position_of_common=common_coeff_to_fork.index(plus_or_minus_one)
            temp_list=deepcopy(list(root_list[-2].coeffs))
            temp_list.pop(position_of_common)
            plus_or_minus_one2=[el for el in temp_list if el!=0][0]
            if plus_or_minus_one*plus_or_minus_one2>0:
                #node -2 is S- node, invert it with S+
                temp_root=root_list_cp[-2]
                root_list_cp[-2]=root_list_cp[-1]
                root_list_cp[-1]=temp_root
        return root_list_cp

class Dynkin:
    def __init__(self,v_arr,rootlist):
        self.v_arr=v_arr
        self.vect_list=v_arr.vect_list
        self.root_list=reorder(rootlist)
        self.roots_coeffs=self.roots_coeffs()
        self.Qindices = self.Q_indices_reorder()
        self.Qvector=self.Q_vector_reps()
        self.reps=self.reps()
        self.cartan_matrix=self.cartan_matrix()
    
    def so_type(self):
        if 2 not in self.root_list[-1].coeffs and -2 not in self.root_list[-1].coeffs:
            return True
        else:
            return False
    
    def write_indices_Q(self,list_indices,root_number):
        #no long roots
        if self.so_type():
            string='Q(' if root_number<sum(self.v_arr.bosons_fermions)-2 else 'S('
            for index in list_indices[0]:
                string+=index
            string+=')|('
            for index in list_indices[1]:
                string+=index
            return string+')'
        else:
            string='Q(' if root_number<sum(self.v_arr.bosons_fermions)-1 else 'S('
            if root_number!=sum(self.v_arr.bosons_fermions)-2:
                for index in list_indices[0]:
                    string+=index
                string+=')|('
                for index in list_indices[1]:
                    string+=index
            else:
                list_indices_sp=self.Qindices[-2]
                list_indices_sm=self.Qindices[-1]
                common_index_bosons=[element for element in list_indices_sp[0] if element in list_indices_sm[0]]
                common_index_fermions=[element for element in list_indices_sp[1] if element in list_indices_sm[1]]
                for index in common_index_bosons:
                    string+=index
                string+=')|('
                for index in common_index_fermions:
                    string+=index
            return string+')'

    def Q_vector_reps(self):
        list=[]
        for i,node in enumerate(self.Qindices):
            list.append(self.write_indices_Q(node,i))
        return list


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
        new_varr=self.v_arr
        new_root_list=[]
        if root_i>len(self.root_list):
            raise ValueError
        refl_root=self.root_list[root_i-1]
        for element in self.root_list:
            new_root_list.append(weyl(refl_root,element))
        return Dynkin(new_varr,reorder(new_root_list))

    def Q_indices(self):
        roots_list=[np.array(root) for root in self.roots_coeffs]
        if self.so_type():
            sum_root=sum(roots_list)-(roots_list[-2]+roots_list[-1])/2
        else:
            sum_root=sum(roots_list)-roots_list[-1]/2
        unique_vectors,indices=[],[]
        bosons=self.v_arr.bosons_fermions[0]
        fermions=self.v_arr.bosons_fermions[1]
        for i in range(bosons+fermions-2):
            unique_vectors.append(list(sum_root))
            sum_root-=roots_list[i]
            indices.append(get_indices_from_unique_vectors(bosons,fermions,indices,unique_vectors[-1]))
        if self.so_type():
            unique_vectors.append(list((roots_list[-2]+roots_list[-1])/2))
            temp_indices=[get_indices_from_unique_vectors(bosons,fermions,indices,unique_vectors[-1])]
            unique_vectors.append(list(-(roots_list[-2]-roots_list[-1])/2))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,unique_vectors[-1]))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,temp_indices,[-el for el in unique_vectors[-1]]))
        else:
            unique_vectors.append(list(sum_root))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,indices,unique_vectors[-1]))
            unique_vectors.append(list(roots_list[-1]/2))
            indices.append(get_indices_from_unique_vectors(bosons,fermions,indices,unique_vectors[-1]))
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

    def QQ_relation(self, node_number):
        node_position=node_number-1
        new_diagram=self.weyl_on_diagram(node_number)
        symm_cartan_row=[]
        for root in self.root_list:
            symm_cartan_row.append(scalarp_roots(self.root_list[node_position],root)/2)
        conn_Qfunctions=[]
        #add Qfunction to RHS only if it's connected by links in Dynkin diagram
        for root_no in range(len(self.root_list)):
            if symm_cartan_row[root_no]!=0 and root_no!=node_position:
                conn_Qfunctions.append(self.Qvector[root_no])
        #LHS is more case by case. GL-type tail is always the same, but for spinor nodes we have to be careful
        if self.so_type(): #so-type diagram
            if node_number==sum(self.v_arr.bosons_fermions)-1: #spinor node +
                if self.root_list[node_position].length!=0: #bosonic QQ
                    product_Q=conn_Qfunctions
                    Wronskian=[self.Qvector[-1],new_diagram.Qvector[-1]] #takes Q functions on spinor - nodes
                    shifts=[symm_cartan_row[node_position],-symm_cartan_row[node_position]]
                else: #fermionic QQ
                    product_Q=[self.Qvector[-1],new_diagram.Qvector[-2]] #takes Q functions on spinor - node and short root of sp diagram
                    Wronskian=[]
                    #add Qfunction to RHS only if it's connected by links in Dynkin diagram
                    for root_no in range(len(self.root_list)):
                        if symm_cartan_row[root_no]!=0 and root_no!=node_position:
                            Wronskian.append(new_diagram.Qvector[root_no])
                    shifts=[self.cartan_matrix[node_position][node_position-1],self.cartan_matrix[node_position][node_position+1]]
            elif node_number==sum(self.v_arr.bosons_fermions): #spinor node -
                if self.root_list[node_position].length!=0: #bosonic QQ
                    product_Q=conn_Qfunctions
                    Wronskian=[self.Qvector[-2],new_diagram.Qvector[-2]] #takes Q functions on spinor + nodes
                    shifts=[symm_cartan_row[node_position],-symm_cartan_row[node_position]]
                else: #fermionic QQ
                    product_Q=[self.Qvector[-2],new_diagram.Qvector[-2]] #takes Q functions on spinor + node and short root of sp diagram
                    Wronskian=[new_diagram.Qvector[-3],new_diagram.Qvector[-1]]
                    shifts=[self.cartan_matrix[node_position][node_position-2],self.cartan_matrix[node_position][node_position-1]]
            else: #gl tail and fork
                if self.root_list[node_position].length!=0: #bosonic QQ
                    product_Q=conn_Qfunctions
                    Wronskian=[self.Qvector[node_position],new_diagram.Qvector[node_position]] 
                    shifts=[2*symm_cartan_row[node_position],-2*symm_cartan_row[node_position]]
                else:
                    Wronskian=conn_Qfunctions
                    product_Q=[self.Qvector[node_position],new_diagram.Qvector[node_position]] 
                    shifts=[self.cartan_matrix[node_position][node_position-1],self.cartan_matrix[node_position][node_position+1]]
        else: #sp-type diagram
            if node_number==sum(self.v_arr.bosons_fermions)-1 and self.root_list[node_position].length==0: #fermionic QQ at short root
                if self.Qvector[-1]==new_diagram.Qvector[-2]:
                    Wronskian=[self.Qvector[-2],new_diagram.Qvector[-1]]
                    product_Q=conn_Qfunctions
                    shifts=[self.cartan_matrix[node_position][node_position-1],self.cartan_matrix[node_position][node_position+1]]
                elif self.Qvector[-1]==new_diagram.Qvector[-1]:
                    Wronskian=[self.Qvector[-2],new_diagram.Qvector[-2]]
                    product_Q=conn_Qfunctions
                    shifts=[self.cartan_matrix[node_position][node_position-1],self.cartan_matrix[node_position][node_position+1]]
            else: #gl tail and fork
                if self.root_list[node_position].length!=0: #bosonic QQ
                    product_Q=conn_Qfunctions
                    Wronskian=[self.Qvector[node_position],new_diagram.Qvector[node_position]] 
                    shifts=[2*symm_cartan_row[node_position],-2*symm_cartan_row[node_position]]
                else:
                    Wronskian=conn_Qfunctions
                    product_Q=[self.Qvector[node_position],new_diagram.Qvector[node_position]] 
                    shifts=[self.cartan_matrix[node_position][node_position-1],self.cartan_matrix[node_position][node_position+1]]
        return {"Wronskian":Wronskian, 'shifts':shifts, 'LHS':product_Q}
    
    def write_QQ(self,node_number):
        QQ_dictionary=self.QQ_relation(node_number)
        Wronskian=QQ_dictionary['Wronskian']
        shifts=QQ_dictionary['shifts']
        if len(Wronskian)==3:
            shifts.append(shifts[-1])
        product_Q=sorted(QQ_dictionary['LHS'])
        LHS='*'.join(product_Q)
        Wronskian_shifts=sorted([Wronskian[i]+f'**[{str(int(shifts[i]))}]' for i in range(len(Wronskian))])
        RHS=','.join(Wronskian_shifts)
        return f'{LHS}=Wronskian[{RHS}]'
    
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

