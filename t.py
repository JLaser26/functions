#Work in progress!

class function:
    def __init__(self, *vars):
        """Add each variable as a string, for ex:- '-2x' """
        self.vars = vars
    
    def __str__(self):
        return " ".join(self.vars)
    
    def tokens(self): return list(self.vars)

    #try to break the function into small parts
    def evaluater(self):
        tk = self.tokens()
        evaluated_tokens = []
        for i in tk:
            for j in i:
                try:
                    x = eval(j)
                except Exception:
                    x = j
                evaluated_tokens.append(x)
        return evaluated_tokens
    
    def functionality(self):
        ...
    
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
y = function("x", "- 2w", "+ 3z")
print(y, y.number_of_variables())
print(y.evaluater())