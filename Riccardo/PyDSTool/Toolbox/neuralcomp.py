#Author:    Robert Clewley
#Date:      09 May 2006
#$Revision: 1.2 $
"""
    An example set of basic compartmental model ModelSpec classes
    for use in computational neuroscience modeling.

    Basic classes provided:
    channel, compartment

    For multi-compartment models:
    soma, dendr_compartment, synapse, neuron, network, exc_synapse, inh_synapse

Rob Clewley, September 2005.
"""

# version 1.2 -- added support for facilitating synapses

from PyDSTool import *
from copy import copy

# which generators will specs made from these classes be compatible with
compatGens = findGenSubClasses('ODEsystem')

# -----------------------------------------------------------------------------

# basic classes -- not intended for direct use
class channel(LeafComponent):
    pass

class compartment(Component):
    pass


# use these for multi-compartment cell model:
class soma(compartment):
    pass

class dendr_compartment(compartment):
    pass

class neuron(Component):
    pass

class synapse(LeafComponent):
    pass

class exc_synapse(synapse):
    pass

class inh_synapse(synapse):
    pass

class network(Component):
    pass

global voltage, V
voltage = 'V'
V = Var(voltage)

# ----------------------------------------------------------------------------
## Object factory functions, and related utilities

# Synapse (dis-)connection

def disconnectSynapse(syn, mspec):
    """Disconnect synapse object syn from ModelSpec object mspec in
    which it has been declared, and remove its target cell's synaptic
    channels."""
    targsfound = []
    for targname in syn.connxnTargets:
        targsfound.extend(searchModelSpec(mspec, targname))
    if len(targsfound) == len(syn.connxnTargets):
        if remain(targsfound, syn.connxnTargets) == []:
            # all targets were found
            mspec.remove(targsfound)
        else:
            raise ValueError("synapse targets were not found in model spec"
                             "argument")
    else:
        raise ValueError("synapse targets were not found in model spec"
                         "argument")
    syn.delConnxnTarget(targsfound)
    mspec.remove(syn)


def connectWithSynapse(synname, syntypestr, source_cell, dest_cell,
                       dest_compartment_name="",
                       threshfun=None, alpha=None, beta=None,
                       threshfun_d=None, alpha_d=None, beta_d=None,
                       adapt_typestr=None,
                       vrev=None, g=None, noauxs=True):
    """Valid source_cell and dest_cell is a soma or a neuron.
    The optional argument dest_compartment_name is a declared name in
    dest_cell, or by default will be the soma.
    """
    if isinstance(source_cell, soma):
        addedNeuron1Name = ""
        soma1name = source_cell.name
    elif isinstance(source_cell, neuron):
        addedNeuron1Name = source_cell.name + '.'
        soma1name = source_cell._componentTypeMap['soma'][0].name
    else:
        raise TypeError("Invalid cell type to connect with synapse")
    if isinstance(dest_cell, soma):
        addedNeuron2Name = ""
        if dest_compartment_name == "":
            comp2name = dest_cell.name
        else:
            raise ValueError("Cannot specify destination compartment name "
                             "when dest_cell has type soma")
    elif isinstance(dest_cell, neuron):
        addedNeuron2Name = dest_cell.name + '.'
        if dest_compartment_name == "":
            comp2name = dest_cell._componentTypeMap['soma'][0].name
        else:
            comp_objlist = dest_cell._componentTypeMap['compartment']
            comp_names = [o.name for o in comp_objlist]
            # check legitimacy of dest_compartment_name in dest_cell
            if dest_cell._componentTypeMap['soma'][0].name == dest_compartment_name:
                comp2name = dest_compartment_name
            elif dest_compartment_name in comp_names:
                comp2name = dest_compartment_name
            else:
                raise ValueError("Destination compartment name " \
                                 + dest_compartment_name \
                                 + "not found in destination cell")
    else:
        raise TypeError("Invalid cell type to connect with synapse")
    if addedNeuron1Name != "":
        # we were given a neuron for source_cell, so should use its name
        # to help distinguish multiple occurrences of svar name at time
        # of registration in a network component.
        svarname1 = source_cell.name+"_"+soma1name
    else:
        svarname1 = soma1name
    if addedNeuron2Name != "":
        # we were given a neuron for dest_cell, so should use its name
        # to help distinguish multiple occurrences of svar name at time
        # of registration in a network component.
        svarname2 = dest_cell.name+"_"+comp2name
    else:
        svarname2 = comp2name
    svar = Var('s_'+svarname1+'_'+svarname2)
    svarname = nameResolver(svar, dest_cell)
    channel_s12 = makeSynapseChannel(synname,
                                     addedNeuron1Name+synname+'.'+svarname,
                                     comp2name+'.'+voltage,
                                     syntypestr, vrev=vrev, g=g, noauxs=noauxs)
    channel_name = nameResolver(channel_s12, dest_cell)
    # update with globalized name
    targetChannel = addedNeuron2Name+comp2name+'.'+channel_name
    channel_s12.name = channel_name
    # check that this source synaptic gating variable hasn't been made before,
    # otherwise re-use that variable?
    if threshfun_d == alpha_d == beta_d == None:
        syn12 = makeSynapse(synname, svarname,
                        addedNeuron1Name+soma1name+'.'+voltage,
                        syntypestr, threshfun=threshfun, alpha=alpha,
                        beta=beta, targetchannel=targetChannel, noauxs=noauxs)
    else:
        syn12 = makeAdaptingSynapse(synname, svarname, 'd',
                        addedNeuron1Name+soma1name+'.'+voltage,
                        syntypestr, adapt_typestr,
                        threshfun=threshfun, alpha=alpha,  beta=beta,
                        threshfun_d=threshfun_d, alpha_d=alpha_d, beta_d=beta_d,
                        targetchannel=targetChannel, noauxs=noauxs)
    if addedNeuron2Name == "":
        # we were given a soma for dest_cell
        dest_cell.add(channel_s12)
    else:
        comp = dest_cell.components[comp2name]
        dest_cell.remove(comp2name)
        comp.add(channel_s12)
        dest_cell.add(comp)
    if addedNeuron1Name != "":
        # we were given a neuron class for source_cell, so we have to
        # add the synaptic variable component to that cell
        source_cell.add(syn12)
    return syn12


# ----------------------------------------------------------------------------


def makeSynapseChannel(name, gatevarname, voltage=voltage, typestr=None,
                       vrev=None, g=None, noauxs=True):
    # make a chemical synapse channel
    channel_s = channel(name, targetLangs, compatibleGens=compatGens,
                              compatibleNodes=[compartment])
    if vrev is None:
        if typestr == 'exc':
            vrevpar = Par('0', 'vrev')
        elif typestr == 'inh':
            vrevpar = Par('-75', 'vrev')
        elif isinstance(typestr, str):
            vrevpar = Par('vrev')
        else:
            raise TypeError("Invalid type for synapse")
    elif isinstance(vrev, str):
        vrevpar = Par(vrev, 'vrev')
    elif isinstance(vrev, float) or isinstance(vrev, int):
        vrevpar = Par(repr(vrev), 'vrev')
    else:
        # allow for ModelSpec function of other variables
        vrevpar = vrev
    v_pot = Var(voltage)   # this is to provide the local name only
    if g is None:
        gpar = Par('g')
    elif isinstance(g, str):
        gpar = Par(g, 'g')
    else:
        gpar = Par(repr(g), 'g')
    gatevar = Var(gatevarname)
    condQ = gpar*gatevar

    if noauxs:
        I = Var(condQ*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        channel_s.add([gpar,vrevpar,I])
    else:
        cond = Var(condQ, 'cond', [0,Inf], specType='ExpFuncSpec')
        shunt = Var(cond*vrevpar, 'shunt', specType='ExpFuncSpec')
        I = Var(cond*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        channel_s.add([gpar,vrevpar,I,cond,shunt])
    return channel_s


def makeBiasChannel(name, I=None, noauxs=True):
    # bias current "channel"
    channel_Ib = channel(name, targetLangs, compatibleGens=compatGens,
                              compatibleNodes=[compartment])
    if I is None:
        Ibias = Par('Ibias')
    elif isinstance(I, str):
        Ibias = Par(I, 'Ibias')
    else:
        Ibias = Par(repr(I), 'Ibias')
    Ibias_cond = Var('0', 'cond', [0,Inf])   # no conductance for this plain current
    Ib = Var(-Ibias, 'I')
    Ibias_shunt = Var(Ib, 'shunt', specType='ExpFuncSpec')    # technically not a shunt but acts like a V-independent shunt!
    if noauxs:
        channel_Ib.add([Ibias,Ib])
    else:
        channel_Ib.add([Ibias,Ibias_cond,Ibias_shunt,Ib])
    return channel_Ib


def makeChannel_halfact(name,voltage=voltage,s=None,isinstant=False,sinf=None,taus=None,spow=1,
                s2=None,isinstant2=False,sinf2=None,taus2=None,spow2=1,vrev=None,g=None,
                parlist=[],noauxs=True):
    # a regular ionic membrane channel
    #
    # assumes that arate's, brate's are functions of the voltage, in a
    # differential equation for s that has the form:
    #   s' = (sinf(v) - v)/taus(v)
    # channel may have up to two gating variables, each of which is
    # given by an ODE
    if g is None:
        gpar = Par('g')
    elif isinstance(g, str):
        gpar = Par(g, 'g')
    else:
        gpar = Par(repr(g), 'g')
    if vrev is None:
        vrevpar = Par('vrev')
    elif isinstance(vrev, str):
        vrevpar = Par(vrev, 'vrev')
    elif isinstance(vrev, float) or isinstance(vrev, int):
        vrevpar = Par(repr(vrev), 'vrev')
    else:
        # allow for ModelSpec function of other variables
        vrevpar = vrev
    declare_list = [gpar,vrevpar]
    thischannel = channel(name, targetLangs, compatibleGens=compatGens,
                          compatibleNodes=[compartment])
    v_pot = Var(voltage)   # this is to provide the local name only
    # deal with gating variable, if present
    if s is None:
        condQ = gpar
        assert taus==sinf==None
        assert spow==1, "power should be unused when gating variable is None"
        # override
        spow=0
    else:
        if taus==sinf==None:
            assert isinstant == False
            gatevar = Par(s)
        else:
            if isinstance(sinf, QuantSpec):
                sinf_term = Var(sinf(), name=s+'inf')
            elif isinstance(sinf, str):
                sinf_term = Var(sinf, name=s+'inf')
            else:
                raise TypeError("sinf argument must be a string or a QuantSpec")
            if isinstant:
                taus_term = Var('0', name='tau'+s)
                if noauxs:
                    gatevar = Var(sinf_term(), name=s, domain=[0,1],
                              specType='ExpFuncSpec')
                else:
                    # re-use aterm and bterm defs for efficiency
                    gatevar = Var(sinf_term(), name=s, domain=[0,1],
                              specType='ExpFuncSpec')
            else:
                if isinstance(taus, QuantSpec):
                    taus_term = Var(taus(), name='tau'+s)
                elif isinstance(taus, str):
                    taus_term = Var(taus, name='tau'+s)
                else:
                    raise TypeError("taus argument must be a string or a QuantSpec")
                # temporary declaration of symbol s (the argument string) as a Var
                s_ = Var(s)
                if noauxs:
                    gatevar = Var( (sinf_term() - s_) / taus_term(), name=s,
                              domain=[0,1],
                              specType='RHSfuncSpec')
                else:
                    # re-use aterm and bterm defs for efficiency
                    gatevar = Var( (sinf_term() - s_) / taus_term(), name=s,
                              domain=[0,1],
                              specType='RHSfuncSpec')
            if not noauxs:
                declare_list.extend([taus_term,sinf_term])
        declare_list.append(gatevar)
        condQ = gpar*makePowerSpec(gatevar,spow)

    # deal with second gating variable, if present
    if s2 is not None:
        assert s is not None, "Must use first gating variable name first!"
        assert spow2 != 0, "Do not use second gating variable with spow2 = 0!"
        s2_ = Var(s2)
        if isinstance(sinf2, QuantSpec):
            sinf2_term = Var(sinf2(), name=s2+'inf', specType='ExpFuncSpec')
        elif isinstance(sinf2, str):
            sinf2_term = Var(sinf2, name=s2+'inf', specType='ExpFuncSpec')
        else:
            raise TypeError("sinf2 argument must be a string or a QuantSpec")
        if isinstant2:
            taus2_term = Var('0', name='tau'+s2)
            if noauxs:
                gatevar2 = Var(sinf2_term(), name=s2, domain=[0,1],
                           specType='ExpFuncSpec')
            else:
                gatevar2 = Var(sinf2_term(), name=s2, domain=[0,1],
                           specType='ExpFuncSpec')
        else:
            if isinstance(taus2, QuantSpec):
                taus2_term = Var(taus2(), name='tau'+s2, specType='ExpFuncSpec')
            elif isinstance(taus2, str):
                taus2_term = Var(taus2, name='tau'+s2, specType='ExpFuncSpec')
            else:
                raise TypeError("taus2 argument must be a string or a QuantSpec")
            if noauxs:
                gatevar2 = Var( (sinf2_term() - s2_) / taus2_term(), name=s2,
                           domain=[0,1], specType='RHSfuncSpec')
            else:
                gatevar2 = Var( (sinf2_term - s2_) / taus2_term(), name=s2,
                           domain=[0,1], specType='RHSfuncSpec')
        declare_list.append(gatevar2)
        if not noauxs:
            declare_list.extend([taus2_term, sinf2_term])
        condQ = condQ * makePowerSpec(gatevar2, spow2)

    if noauxs:
        I = Var(condQ*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        declare_list.append(I)
    else:
        cond = Var(condQ, 'cond', [0,Inf], specType='ExpFuncSpec')
        shunt = Var(cond*vrevpar, 'shunt', specType='ExpFuncSpec')
        I = Var(cond*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        declare_list.extend([I,cond,shunt])
    declare_list.extend(parlist)
    thischannel.add(declare_list)
    return thischannel


def makeChannel_rates(name,voltage=voltage,
                      s=None,isinstant=False,arate=None,brate=None,spow=1,
                      s2=None,isinstant2=False,arate2=None,brate2=None,spow2=1,
                      vrev=None,g=None,parlist=[],noauxs=True):
    # a regular ionic membrane channel
    #
    # assumes that arate's, brate's are functions of the voltage, in a
    # differential equation for s that has the form:
    #   s' = sa(v) * (1-s) - sb(v) * s
    # channel may have up to two gating variables, each of which is
    # given by an ODE
    if g is None:
        gpar = Par('g')
    elif isinstance(g, str):
        gpar = Par(g, 'g')
    else:
        gpar = Par(repr(g), 'g')
    if vrev is None:
        vrevpar = Par('vrev')
    elif isinstance(vrev, str):
        vrevpar = Par(vrev, 'vrev')
    elif isinstance(vrev, float) or isinstance(vrev, int):
        vrevpar = Par(repr(vrev), 'vrev')
    else:
        # allow for ModelSpec function of other variables
        vrevpar = vrev
    declare_list = [gpar,vrevpar]
    thischannel = channel(name, targetLangs, compatibleGens=compatGens,
                          compatibleNodes=[compartment])
    v_pot = Var(voltage)   # this is to provide the local name only
    # deal with gating variable, if present
    if s is None:
        condQ = gpar
        assert arate==brate==None
        assert spow==1, "power should be unused when gating variable is None"
        # override
        spow=0
    else:
        if arate==brate==None:
            assert isinstant == False
            gatevar = Par(s)
        else:
            if isinstance(arate, QuantSpec):
                aterm = Var(arate(), name='a'+s)
            elif isinstance(arate, str):
                aterm = Var(arate, name='a'+s)
            else:
                raise TypeError("arate argument must be a string or QuantSpec")
            if isinstance(brate, QuantSpec):
                bterm = Var(brate(), name='b'+s)
            elif isinstance(brate, str):
                bterm = Var(brate, name='b'+s)
            else:
                raise TypeError("brate argument must be a string or QuantSpec")
            # temporary declaration of symbol s (the argument string) as a Var
            s_ = Var(s)
            if isinstant:
                taus_term = Var('0', name='tau'+s)
                sinf_term = Var( aterm/(aterm+bterm), name=s+'inf')
                if noauxs:
                    gatevar = Var( aterm()/(aterm()+bterm()), name=s,
                              domain=[0,1], specType='ExpFuncSpec')
                else:
                    # re-use aterm and bterm defs for efficiency
                    gatevar = Var( aterm()/(aterm()+bterm()), name=s,
                              domain=[0,1], specType='ExpFuncSpec')
            else:
                taus_term = Var( 1/(aterm+bterm), name='tau'+s)
                sinf_term = Var( aterm*taus_term, name=s+'inf')
                if noauxs:
                    gatevar = Var( aterm()*(1-s_)-bterm()*s_, name=s,
                              domain=[0,1], specType='RHSfuncSpec')
                else:
                    # re-use aterm and bterm defs for efficiency
                    gatevar = Var( aterm()*(1-s_)-bterm()*s_, name=s,
                              domain=[0,1], specType='RHSfuncSpec')
            if not noauxs:
                declare_list.extend([taus_term,sinf_term,aterm,bterm])
        declare_list.append(gatevar)
        condQ = gpar*makePowerSpec(gatevar,spow)

    # deal with second gating variable, if present
    if s2 is not None:
        assert s is not None, "Must use first gating variable name first!"
        assert spow2 != 0, "Do not use second gating variable with spow2 = 0!"
        if isinstance(arate2, QuantSpec):
            aterm2 = Var(arate2(), name='a'+s2)
        elif isinstance(arate2, str):
            aterm2 = Var(arate2, name='a'+s2)
        else:
            raise TypeError("arate2 argument must be a string or QuantSpec")
        if isinstance(brate2, QuantSpec):
            bterm2 = Var(brate2(), name='b'+s2)
        elif isinstance(brate2, str):
            bterm2 = Var(brate2, name='b'+s2)
        else:
            raise TypeError("brate2 argument must be a string or QuantSpec")
        s2_ = Var(s2)
        if isinstant2:
            if noauxs:
                gatevar2 = Var( aterm2()/(aterm2()+bterm2()), name=s2,
                              domain=[0,1], specType='ExpFuncSpec')
            else:
                taus2_term = Var( '0', name='tau'+s2)
                sinf2_term = Var( aterm2()/(aterm2()+bterm2()), name=s2+'inf')
                gatevar2 = Var( aterm2()/(aterm2()+bterm2()), name=s2,
                              domain=[0,1], specType='ExpFuncSpec')
        else:
            if noauxs:
                gatevar2 = Var(aterm2()*(1-s2_)-bterm2()*s2_, name=s2,
                           domain=[0,1], specType='RHSfuncSpec')
            else:
                taus2_term = Var( 1/(aterm2()+bterm2()), name='tau'+s2)
                sinf2_term = Var( aterm2()*taus2_term, name=s2+'inf')
                gatevar2 = Var(aterm2()*(1-s2_)-bterm2()*s2_, name=s2,
                           domain=[0,1], specType='RHSfuncSpec')
        condQ = condQ * makePowerSpec(gatevar2, spow2)
        declare_list.append(gatevar2)
        if not noauxs:
            declare_list.extend([taus2_term, sinf2_term, aterm2, bterm2])

    if noauxs:
        I = Var(condQ*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        declare_list.append(I)
    else:
        cond = Var(condQ, 'cond', [0,Inf], specType='ExpFuncSpec')
        shunt = Var(cond*vrevpar, 'shunt', specType='ExpFuncSpec')
        I = Var(cond*(v_pot-vrevpar), 'I', specType='ExpFuncSpec')
        declare_list.extend([I,cond,shunt])
    declare_list.extend(parlist)
    thischannel.add(declare_list)
    return thischannel


def makeSoma(name, voltage=voltage, channelList=[], C=None, noauxs=True,
             subclass=soma):
    v = Var(QuantSpec(voltage, '-for(channel,I,+)/C', 'RHSfuncSpec'),
        domain=[-130,70])
    if C is None:
        C = Par('C')
    elif isinstance(C, str):
        gpar = Par(C, 'C')
    else:
        C = Par(repr(C), 'C')

    asoma = soma.__new__(subclass)
    asoma.__init__(name, targetLangs, [channel],
                   compatibleNodes=[neuron, network],
                   compatibleGens=compatGens)
    declare_list = copy(channelList)
    if noauxs:
        declare_list.extend([v,C])
    else:
        tauv = Var(QuantSpec('tauv', '1/(C*for(channel,cond,+))'))
        vinf = Var(QuantSpec('vinf', 'for(channel,shunt,+)'+ \
                                       '/for(channel,cond,+)'))
        declare_list.extend([v,C,vinf,tauv])
    asoma.add(declare_list)
    return asoma


def makePointNeuron(name, voltage=voltage, channelList=[], synapseList=[], C=None):
    asoma = makeSoma('soma', voltage, channelList, C=C)
    n = neuron(name, targetLangs, [soma, synapse],
               compatibleNodes=[network], compatibleGens=compatGens)
    n.add(asoma)
    if synapseList != []:
        n.add(synapseList)
    return n


def makePointNeuronNetwork(name, componentList):
    net = network(name, targetLangs, [soma, synapse],
                  compatibleGens=compatGens)
    net.add(componentList)
    return net


def makeNeuronNetwork(name, neuronList):
    # untested!
    net = network(name, targetLangs, [neuron],
                  compatibleGens=compatGens)
    net.add(neuronList)
    return net


def makeSynapse(name, gatevar, precompartment, typestr, threshfun=None,
                alpha=None, beta=None, targetchannel=None,
                evalopt=True, noauxs=True):
    if targetchannel is None:
        raise TypeError("Must provide name of synaptic channel object in "
                        "target cell's compartment")
    if typestr == 'exc':
        if alpha is None:
            alpha = 10.
        if beta is None:
            beta = 0.5
        syn = exc_synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    elif typestr == 'inh':
        if alpha is None:
            alpha = 1.
        if beta is None:
            beta = 0.1
        syn = inh_synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    elif typestr == "":
        syn = synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    else:
        raise ValueError("Invalid type of synapse specified")
    s_ = Var(gatevar)
    if alpha is None:
        a = Par('alpha')
    elif isinstance(alpha, str):
        gpar = Par(alpha, 'alpha')
    else:
        a = Par(repr(alpha), 'alpha')
    if beta is None:
        b = Par('beta')
    elif isinstance(beta, str):
        gpar = Par(beta, 'beta')
    else:
        b = Par(repr(beta), 'beta')
    if threshfun is None:
        f = Fun('0.5+0.5*tanh(v/4.)', ['v'], 'thresh')
    else:
        assert isinstance(threshfun, tuple), \
               "threshfun must be pair (vname, funbody)"
        if isinstance(threshfun[1], QuantSpec):
            funbody = threshfun[1].specStr
        elif isinstance(threshfun[1], str):
            funbody = threshfun[1]
        else:
            raise TypeError("threshold function must be a string or a "
                            "QuantSpec")
        assert threshfun[0] in funbody, \
               "voltage name %s does not appear in function body!"%threshfun[0]
        f = Fun(funbody, [threshfun[0]], 'thresh')
    assert len(f.signature) == 1, \
           'threshold function must be a function of a single argument (voltage)'
    if evalopt:
        alpha_ = a*f.eval(**{f.signature[0]: precompartment})
    else:
        alpha_ = a*f(precompartment)
        syn.add(f)
    s = Var(alpha_*(1-s_)-b*s_, name=gatevar,
            specType='RHSfuncSpec', domain=[0,1])
    syn.add([s,a,b])
    syn.addConnxnTarget(targetchannel)
    if not noauxs:
        sinf = Var( alpha_/(alpha_+b), name=gatevar+'inf')
        taus = Var( 1/(alpha_+b), name='tau'+gatevar)
        syn.add([sinf,taus])
    return syn


def makeAdaptingSynapse(name, gatevar, adaptvar, precompartment, typestr,
                adapt_typestr,
                threshfun=None, alpha=None, beta=None,
                threshfun_d=None, alpha_d=None, beta_d=None,
                targetchannel=None, evalopt=True, noauxs=True):
    if targetchannel is None:
        raise TypeError("Must provide name of synaptic channel object in "
                        "target cell's compartment")
    if typestr == 'exc':
        if alpha is None:
            alpha = 10.
        if beta is None:
            beta = 0.5
        syn = exc_synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    elif typestr == 'inh':
        if alpha is None:
            alpha = 1.
        if beta is None:
            beta = 0.1
        syn = inh_synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    elif typestr == "":
        syn = synapse(name, targetLangs, compatibleNodes=[neuron, network],
                          compatibleGens=compatGens)
    else:
        raise ValueError("Invalid type of synapse specified")
    # synaptic variable
    s_ = Var(gatevar)
    if alpha is None:
        a = Par('alpha')
    elif isinstance(alpha, str):
        a = Par(alpha, 'alpha')
    else:
        a = Par(repr(alpha), 'alpha')
    if beta is None:
        b = Par('beta')
    elif isinstance(beta, str):
        b = Par(beta, 'beta')
    else:
        b = Par(repr(beta), 'beta')
    if threshfun is None:
        f = Fun('0.5+0.5*tanh(v/4.)', ['v'], 'thresh')
    else:
        assert isinstance(threshfun, tuple), \
               "threshfun must be pair (vname, funbody)"
        if isinstance(threshfun[1], QuantSpec):
            funbody = threshfun[1].specStr
        elif isinstance(threshfun[1], str):
            funbody = threshfun[1]
        else:
            raise TypeError("threshold function must be a string or a "
                            "QuantSpec")
        assert threshfun[0] in funbody, \
               "voltage name %s does not appear in function body!"%threshfun[0]
        f = Fun(funbody, [threshfun[0]], 'thresh')
    assert len(f.signature) == 1, \
           'threshold function must be a function of a single argument (voltage)'
    if evalopt:
        alpha_ = a*f.eval(**{f.signature[0]: precompartment})
    else:
        alpha_ = a*f(precompartment)
        syn.add(f)
    # adaptive variable
    if adapt_typestr == 'f':
        at = ''
    elif adapt_typestr == 'd':
        at = '-'
    else:
        raise ValueError, "Invalid type for adapting synapse: use 'f' or 'd'"
    d_ = Var(adaptvar)
    if alpha_d is None:
        a_d = Par('alpha_d')
    elif isinstance(alpha_d, str):
        a_d = Par(alpha_d, 'alpha_d')
    else:
        a_d = Par(repr(alpha_d), 'alpha_d')
    if beta_d is None:
        b_d = Par('beta_d')
    elif isinstance(beta_d, str):
        b_d = Par(beta_d, 'beta_d')
    else:
        b_d = Par(repr(beta_d), 'beta_d')
    if threshfun_d is None:
        f_d = Fun('0.5+0.5*tanh(%sv/4.)'%atype, ['v'], 'thresh')
    else:
        assert isinstance(threshfun_d, tuple), \
               "threshfun must be pair (vname, funbody)"
        if isinstance(threshfun_d[1], QuantSpec):
            funbody_d = threshfun_d[1].specStr
        elif isinstance(threshfun_d[1], str):
            funbody_d = threshfun_d[1]
        else:
            raise TypeError("threshold function must be a string or a "
                            "QuantSpec")
        assert threshfun_d[0] in funbody_d, \
               "voltage name %s does not appear in function body!"%threshfun[0]
        f_d = Fun(funbody_d, [threshfun_d[0]], 'thresh')
    assert len(f_d.signature) == 1, \
           'threshold function must be a function of a single argument (voltage)'
    if evalopt:
        alpha_d_ = a_d*f_d.eval(**{f_d.signature[0]: precompartment})
    else:
        alpha_d_ = a_d*f_d(precompartment)
        syn.add(f_d)
    d = Var(alpha_d_*(1-d_)-b_d*d_, name=adaptvar,
            specType='RHSfuncSpec', domain=[0,1])
    s = Var(alpha_*(d_-s_)-b*s_, name=gatevar,
            specType='RHSfuncSpec', domain=[0,1])
    syn.add([s,a,b,d,a_d,b_d])
    syn.addConnxnTarget(targetchannel)
    if not noauxs:
        sinf = Var( alpha_*d_/(alpha_+b), name=gatevar+'inf')
        taus = Var( 1/(alpha_+b), name='tau'+gatevar)
        dinf = Var( alpha_d_/(alpha_d_+b_d), name=adaptvar+'inf')
        taud = Var( 1/(alpha_d_+b_d), name='tau'+adaptvar)
        syn.add([sinf,taus,dinf,taud])
    return syn

# -----------------------------------------------------------------------------

## Helper functions

def makePowerSpec(v, p):
    if isinstance(p, int):
        if p < 0:
            raise ValueError, "power should be > 0"
        elif p==0:
            return "1"
        elif p > 4:
            return Pow(v, p)
        else:
            res = v
            for i in range(p-1):
                res = res * v
            return res
    else:
        return Pow(v, p)


# the following functions will be removed if the calling order for Quantities
# is changed (as it ought to be some day).

def makePar(parname, val=None):
    if val is None:
        return Par(parname)
    elif isinstance(val, str):
        gpar = Par(val, parname)
    else:
        return Par(repr(val), parname)


def makeFun(funname, sig, defn):
    return Fun(defn, sig, funname)
