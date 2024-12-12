from sympy import *
from sympy.geometry import *

class GeometryDSL:
    def __init__(self):
        self.variables = {}
        self.constraints = []
    
    def define_point(self, name):
        self.variables[name] = Point(Symbol(name + "_x"), Symbol(name + "_y"))
    
    def define_line(self, name, p1, p2):
        self.variables[name] = Segment(self.variables[p1], self.variables[p2])
    
    def define_circle(self, name, center, radius):
        self.variables[name] = Circle(self.variables[center], radius)
    
    def define_angle(self, name, p1, p2, p3):
        p1_pt = self.variables[p1]
        p2_pt = self.variables[p2]
        p3_pt = self.variables[p3]
        v1 = Matrix([p1_pt.x - p2_pt.x, p1_pt.y - p2_pt.y])
        v2 = Matrix([p3_pt.x - p2_pt.x, p3_pt.y - p2_pt.y])
        self.variables[name] = acos((v1.dot(v2)) / (v1.norm() * v2.norm()))
        self.variables[name + "_cos"] = (v1.dot(v2)) / (v1.norm() * v2.norm())
        
    def set_length(self, segment_name, length):
        segment = self.variables[segment_name]
        self.constraints.append(Eq(segment.length, length))
    
    def set_angle(self, angle_name, angle_value):
        self.constraints.append(Eq(self.variables[angle_name + "_cos"], cos(angle_value)))
    
    def add_constraint(self, constraint):
        self.constraints.append(constraint)
    
    def collect_substitutions(self):
        subs_dict = {}
        for eq in self.constraints:
            if isinstance(eq, Eq):
                lhs = eq.lhs
                rhs = eq.rhs
                if lhs.is_Symbol:
                    subs_dict[lhs] = rhs
                elif rhs.is_Symbol:
                    subs_dict[rhs] = lhs
                else:
                    for sym in lhs.free_symbols:
                        subs_dict[sym] = solve(eq, sym)[0]
                    for sym in rhs.free_symbols:
                        subs_dict[sym] = solve(eq, sym)[0]
        return subs_dict
    
    def solve(self, target, assume_positions=True):
        symbols = set()
        for var in self.variables.values():
            if hasattr(var, 'free_symbols'):
                symbols.update(var.free_symbols)
        if assume_positions:
            points = [v for v in self.variables.values() if isinstance(v, Point)]
            fixed_coords = set()
            for c in self.constraints:
                if isinstance(c, Eq):
                    fixed_coords.add(c.lhs)
                    fixed_coords.add(c.rhs)
            if points:
                p0 = points[0]
                if p0.x not in fixed_coords:
                    self.constraints.append(Eq(p0.x, 0))
                if p0.y not in fixed_coords:
                    self.constraints.append(Eq(p0.y, 0))
                if len(points) > 1:
                    p1 = points[1]
                    if p1.y not in fixed_coords:
                        self.constraints.append(Eq(p1.y, 0))
                    if p1.x not in fixed_coords and p1.x != p0.x:
                        self.constraints.append(Ne(p1.x, p0.x)) 
    
        equations = [eq for eq in self.constraints if isinstance(eq, Eq)]
        try:
            solution = solve(equations, symbols, dict=True)
            if solution:
                # Return all solutions
                solutions = solution
                if isinstance(target, str):
                    target_var = self.variables.get(target)
                    if target_var:
                        if isinstance(target_var, Segment):
                            lengths = [target_var.length.subs(sol).evalf() for sol in solutions]
                            return lengths
                        elif not isinstance(target_var, (Point, Circle, Segment)):
                            values = [target_var.subs(sol).evalf() for sol in solutions]
                            return values
                        elif isinstance(target_var, Point):
                            coordinates = []
                            for sol in solutions:
                                x_val = target_var.x.subs(sol).evalf()
                                y_val = target_var.y.subs(sol).evalf()
                                coordinates.append( (x_val, y_val) )
                            return coordinates
                        elif isinstance(target_var, Circle):
                            return [target_var.subs(sol) for sol in solutions]
                        else:
                            values = [target_var.subs(sol).evalf() for sol in solutions]
                            return values
                    else:
                        return "Target variable not found"
                else:
                    # target is an expression
                    values = [target.subs(sol).evalf() for sol in solutions]
                    return values
            else:
                return "No solution found."
        except Exception as e:
            print(f"Error during solving: {e}")
            return "Error during solving."

def example1():
    dsl1 = GeometryDSL()
    dsl1.define_point("A")
    dsl1.define_point("B")
    dsl1.define_point("C")
    dsl1.define_line("AB", "A", "B")
    dsl1.define_line("BC", "B", "C")
    dsl1.define_line("AC", "A", "C")
    dsl1.define_angle("BAC", "B", "A", "C")

    dsl1.set_length("AB", 5)
    dsl1.set_length("BC", 6)
    dsl1.set_length("AC", 7)

    # Set coordinates of points A and B
    dsl1.constraints.append(Eq(dsl1.variables["A"].x, 0))
    dsl1.constraints.append(Eq(dsl1.variables["A"].y, 0))
    dsl1.constraints.append(Eq(dsl1.variables["B"].x, 5))
    dsl1.constraints.append(Eq(dsl1.variables["B"].y, 0))

    # Solve for coordinates of C
    C_solutions = dsl1.solve("C")
    print(f"Example 1: Coordinates of C: {C_solutions}")

    # Let's pick the solution with positive y-coordinate
    C_chosen = None
    for coords in C_solutions:
        C_x, C_y = coords
        if C_y > 0:
            C_chosen = coords
            break
    if C_chosen is None:
        # If none have positive y-coordinate, use any solution
        C_chosen = C_solutions[0]

    C_x, C_y = C_chosen
    print(f"Using Coordinates of C: ({C_x}, {C_y})")

    # Compute angle BAC using Law of Cosines
    A_x, A_y = 0, 0
    B_x, B_y = 5, 0

    # Vectors BA and CA
    BA = Matrix([A_x - B_x, A_y - B_y])
    CA = Matrix([A_x - C_x, A_y - C_y])

    # Compute angle at A between vectors BA and CA
    cos_angle_BAC = (BA.dot(CA)) / (BA.norm() * CA.norm())
    angle_BAC = acos(cos_angle_BAC)
    print(f"Example 1: Angle BAC: {angle_BAC.evalf()}")

def example2():
    dsl2 = GeometryDSL()
    dsl2.define_point("O")
    dsl2.define_point("P")
    dsl2.define_circle("c", "O", 5)

    dsl2.constraints.append(Eq(dsl2.variables["O"].x, 0))
    dsl2.constraints.append(Eq(dsl2.variables["O"].y, 0))
    dsl2.constraints.append(Eq(dsl2.variables["P"].x, 3))
    dsl2.constraints.append(Eq(dsl2.variables["P"].y, 4))

    subs_dict = dsl2.collect_substitutions()

    circle = dsl2.variables["c"]
    point = dsl2.variables["P"]

    distance_expr = sqrt((point.x - circle.center.x)**2 + (point.y - circle.center.y)**2)
    distance_value = distance_expr.subs(subs_dict).evalf()
    O_x_val = circle.center.x.subs(subs_dict).evalf()
    O_y_val = circle.center.y.subs(subs_dict).evalf()
    P_x_val = point.x.subs(subs_dict).evalf()
    P_y_val = point.y.subs(subs_dict).evalf()

    print(f"Example 2: Circle equation: Center({O_x_val}, {O_y_val}), Radius: {circle.radius}")
    print(f"Example 2: Distance from O to P: {distance_value}")
    print(f"Example 2: Point P inside circle: {distance_value <= circle.radius}")

def example3():
    dsl3 = GeometryDSL()
    dsl3.define_point("A")
    dsl3.define_point("B")
    dsl3.define_point("C")
    dsl3.define_line("AB", "A", "B")
    dsl3.define_line("BC", "B", "C")
    dsl3.define_line("AC", "A", "C")
    dsl3.define_angle("BAC", "B", "A", "C")
    dsl3.define_angle("ABC", "A", "B", "C")
    dsl3.define_angle("BCA", "B", "C", "A")

    dsl3.set_length("AB", 5)
    dsl3.set_length("BC", 5)
    dsl3.set_length("AC", 5)

    # Set coordinates of points A and B
    dsl3.constraints.append(Eq(dsl3.variables["A"].x, 0))
    dsl3.constraints.append(Eq(dsl3.variables["A"].y, 0))
    dsl3.constraints.append(Eq(dsl3.variables["B"].x, 5))
    dsl3.constraints.append(Eq(dsl3.variables["B"].y, 0))

    C_solutions = dsl3.solve("C")
    print(f"Example 3: Coordinates of C: {C_solutions}")


if __name__ == '__main__':
    # LLM generated examples:
    example1()
    example2()
    example3()