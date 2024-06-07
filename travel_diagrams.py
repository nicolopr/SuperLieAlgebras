def weyl_transf(diagram,node):
    diagram2=diagram.weyl_on_diagram(node)
    return diagram2

def travel(diagram):
    visited_diagrams, roots_reps=[],[]
    diagrams_to_transform=[diagram]
    while diagrams_to_transform:
        to_transform=diagrams_to_transform.pop()
        for i in range(1,sum(diagram.v_arr.bosons_fermions)+1):
            new_diag=weyl_transf(to_transform,i)
            if new_diag.roots_coeffs not in roots_reps:
                visited_diagrams.append(new_diag)
                roots_reps.append(new_diag.roots_coeffs)
                diagrams_to_transform.append(new_diag)
    return visited_diagrams

def find_Qfunctions(visited_diagrams):
    Q_functions=[]
    for diagram in visited_diagrams:
        skip=False
        if 2 in diagram.root_list[-1].coeffs or -2 in diagram.root_list[-1].coeffs:
            #skip long roots diagram
            skip=True
        for Q_function in diagram.Qvector:
            if Q_function not in Q_functions and not skip:
                Q_functions.append(Q_function)
    return Q_functions 

def find_QQ_relations(visited_diagrams):
    QQ_rels=[]
    for diagram in visited_diagrams:
        for node_no in range(1,4):
            if diagram.QQ_relation(node_no) not in QQ_rels:
                QQ_rels.append(diagram.QQ_relation(node_no))
    return QQ_rels

def print_visited(visited):
    for element in visited:
        element.view()

def print_visited_cartan(visited):
    for element in visited:
        element.view()
        element.print_cartan()