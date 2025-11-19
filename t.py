#Work in progress!

class function:
    def __init__(self, *vars):
        """Add each variable as a string, for ex:- '-2x' """
        self.vars = vars
    
    def __str__(self):
        return " ".join(self.vars)
    
    def functionality(self):
        return
    
    def output(self): return self.functionality()

    def Domain(self):
        ...

    def Range(self):
        ...
    
    def number_of_variables(self): return len(self.vars)

    def derivative(self):
        ...
    
    def integral(self):
        ...

#testing
y = function("x", "- 2y", "+ 3z")
print(y, y.number_of_variables())