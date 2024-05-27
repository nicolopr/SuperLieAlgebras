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

def print_visited(visited):
    for element in visited:
        element.view()

def print_visited_cartan(visited):
    for element in visited:
        element.view()
        element.print_cartan()