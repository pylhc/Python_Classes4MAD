from math import log10

pre={}
pre[0] =''
pre[1] ='k'
pre[2] ='M'
pre[3] ='G'
pre[4] ='T'
pre[5] ='P'
pre[6] ='E'
pre[-1] ='m'
pre[-2] ='u'
pre[-3] ='n'
pre[-4] ='p'
pre[-5] ='f'
pre[-6] ='a'

def num_to_str (n,d=3,p=True,tex=0):
    '''Convert a number to a string in engineering notation.  E.g., 5e-9 -> 5n'''
    if tex:
      pre[-2] ='\\mu'
    if n==0:
      man='0'
      suffix=''
    else:
      order =int(log10(abs(n)))
      if (type(n) is int) and (order<3):
        man="%i" % n
        suffix=''
      else:
        man='%.'+str(2-order%3)+'f'
        man=man % (n / 10.**(3*(order/3)))
        suffix='e'+str(3*(order/3))
        if p and abs(order/3)<7:
          suffix=pre[order/3]
    if tex:
      return (man,suffix)
    else:
      return man+suffix
