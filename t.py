class function:
    def __init__(self, *vars):
        self.vars = vars
    
    def functionality(self, vars):
        return
    
    def output(self): return self.functionality()

    def Domain(self):
        ...

    def Range(self):
        ...
    
    def number_of_variables(self): return len(self.vars)

    def derivative(self) -> function:
        ...
    
    def integral(self) -> function:
        ...

