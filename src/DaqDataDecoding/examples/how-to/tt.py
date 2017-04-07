import math
import ROOT

def get_positions(nt,n_trigs=16,cut=0.2,len_min=10):
    
    vals={}
    for i in range(n_trigs):
        vals[i] = []
    
    for entry in range(nt.GetEntries()):
        nt.GetEntry(entry)
        hits = nt.GetArgs()[1]
        mt = nt.GetArgs()[2]
        if hits!=1:
            continue    # Analyse only pure triggers!
        tt=[]
        for i in range(n_trigs):
            v = nt.GetArgs()[3+i]
            if v!=0:
                vals[i].append(v)

    pos={}
    sigma={}
    for bit,array in vals.iteritems():
        pos[bit]=0
        sigma[bit]=0
        array.sort()
        if len(array)<len_min:
            print 'too small entries (%d) for trigger bit %d!' % (len(array),bit)
            continue
        # print bit,len(array)

        # Cut 20% of entries from begin and end.
        n = int(len(array)*cut)
        array_cut = array[n:-n]
        assert (len(array_cut)+2*n)==len(array)
        m = len(array_cut)
        pos[bit] = sum(array_cut)/m
        sum_sqr = sum([q*q for q in array_cut])
        
        sigma[bit] = math.sqrt(sum_sqr/(m-1)-(pos[bit])**2)
        
        print 'bit:%2d   mean = %g   sigma = %g' % (bit,pos[bit],sigma[bit])
    
    return pos,sigma

def draw_tt(tt,pos,sigma):

    ROOT.gStyle.SetOptFit(1111)
    hists = []
    c= ROOT.TCanvas('tt','tt')
    c.Divide(4,4)
    for i in range(16):
        c.cd(i+1)
        what = 'tt%d' % i
        
        h = ROOT.TH1F('h_tt%d' % i,'h',100,pos[i]-10,pos[i]+10)
        hists.append(h)
        tt.Draw(what+'>>'+h.GetName(),'hits==1&&abs(%s-%g)<10' % (what,pos[i]) )
        h.Fit('gaus')

    c.cd(0)
    c.Print()

def analyze(nt,pos,n_trigs=16):
    for entry in range(nt.GetEntries()):
        nt.GetEntry(entry)
        hits = nt.GetArgs()[1]
        mt = nt.GetArgs()[2]

        tt_corr={}
        for i in range(n_trigs):
            v = nt.GetArgs()[3+i]
            if v:
                tt_corr[v-pos[i]] = i
        if len(tt_corr)>1:
            tt_sort = tt_corr.keys()
            tt_sort.sort()
            for t in tt_sort:
                bit = tt_corr[t]
                print '(%2d,%+5.1f) ' % (bit,t),
            print
        
if __name__=='__main__':

    f = ROOT.TFile('tt2.root')
    tt = f.Get('tt')

    pos,sigma = get_positions(tt)

    draw_tt(tt,pos,sigma)

    #analyze(tt,pos)
