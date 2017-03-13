#!/usr/bin/env python
## from https://opensource.googleblog.com/2017/03/python-fire-command-line.html?m=1

import fire
import scipy
import scipy.optimize

class Example(object):
    def hello_fellas(self, names=['earth','mars','venus']):
        """ Say hello to a list of specified names."""
        return ['Hello {name}!'.format(name=name) for name in names]

def main():
    fire.Fire(Example)

if __name__ == '__main__':
    main()
    
#class Solver(object):
#    def nonlinear(self):#, eq=[x[0]**2 - 2*x[0]*x[1] - x[2], x[0]*x[1] + x[1]**2*x[2], x[2]*x[0] - x[1]*x[2]+ x[0] - 1]):
#        """Solve a system of non-linear equations."""
#        f = lambda x:eq
#        x0 = scipy.optimize.fsolve(f, [1, -1, 2])
#        return x0


#def main():
#    fire.Fire(Solver)

#if __name__ == '__main__':
#    main()



#>>> import scipy
#>>> import scipy.optimize
#>>> f = lambda x: [x[0]**2 - 2*x[0]*x[1] - x[2], x[0]*x[1] + x[1]**2*x[2], x[2]*x[0] - x[1]*x[2]+ x[0] - 1]    
#>>> x0 = scipy.optimize.fsolve(f, [1, -1, 2])
#>>> print x0
