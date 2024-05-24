from roots import *
from weyl import weyl
from basis_vectors import *

class Dynkin:
    def __init__(self,v_arr,root_list):
        self.v_arr=v_arr
        self.vect_list=v_arr.vect_list
        self.root_list=root_list
        self.reps=self.reps()
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
    def view(self):
        print(self.reps[0])
        print(self.reps[1])

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
v_arr=Vector_array((2,1))
print(v_arr)
dk=Distinguished(v_arr)
dk.view()