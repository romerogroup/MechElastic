from mechelastic import EOS

eos = EOS("PvsV.dat", eostype="pressure")
eos.plot_eos()
