from ROOT import TF1, TCanvas

def identity( x ):
    return x[0]

# create an identity function
f = TF1('pyf1', identity, -1., 1.)

# plot the function
c = TCanvas()
f.Draw()
