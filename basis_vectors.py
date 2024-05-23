class Basis_vector:
    '''
    Basic class. Stores grading and index of a basis vector of the root space.
    '''
    def __init__(self,index: int,grading: int):
        self.grading=grading
        self.index=index
    def reps(self):
        if self.grading==0:
            return f'b{self.index}'
        else:
            return f'f{self.index}'
    def view(self):
        print(self.reps())

def vect_list(bosons_fermions: tuple[int,int]):
        '''
    Auxiliary function returning a list of basis vector objects given a tuple number of bosonic and fermionic directions. Used in Vector_array class
    '''
        vect_list=[]
        for i in range(1,bosons_fermions[0]+1):
            vect_list.append(Basis_vector(i,0))
        for i in range(1,bosons_fermions[1]+1):
            vect_list.append(Basis_vector(i,1))
        return vect_list

class Vector_array:
    '''
    Creates a list of basis vectors of root space given given a tuple number of bosonic and fermionic directions.
    '''
    def __init__(self, bosons_fermions: tuple[int,int]):
        self.bosons_fermions=bosons_fermions
        self.vect_list=vect_list(bosons_fermions)
    def reps(self):
        res=''
        for vect in self.vect_list:
            res+=vect.reps()+', '
        res.strip(', ')
        return res
    def view(self):
        print(self.reps())