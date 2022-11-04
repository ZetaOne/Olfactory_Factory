TITLE neurorishika_excitatory.mod   Kenyon cells from wustenberg 2004
 
COMMENT

 : Filename: nas_wustenberg.mod
: Description: 
: Author: Subhasis Ray
: Maintainer: 
: Created: Wed Dec 13 19:06:03 EST 2017
: Version: 
: Last-Updated: Mon Jun 18 14:38:15 2018 (-0400)
:           By: Subhasis Ray
: URL: 
: Doc URL: 
: Keywords: 
: Compatibility: 
: 
: 

: Commentary: 
: 
: NEURON implementation of slow Na+ channel ( NAS ) from Wustenberg
: DG, Boytcheva M, Grunewald B, Byrne JH, Menzel R, Baxter DA


TITLE Honey bee KC from Wustenberg et al 2004
 ***************************************************************************

 changes :
 included the leak current in KCs 

  - USEION removed as we don't want to deal with ions and set eNa and
  eK directly. Rev potentials 'egna' and 'egk' are in the PARAMETERS
  list
    
  - temp: set to 6.3 Celsius, alpha and beta values are set/manipulated
  directly to simulate characteristic firing pattern

  I. Vida, Nov. 2000

  ***************************************************************************
ENDCOMMENT
 
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

? interface

NEURON {
        SUFFIX hh_wuKC
        NONSPECIFIC_CURRENT inaf,inas,ikv,ia,il,ist

        RANGE gnafbar,gnaf, egnaf , gnasbar,gnas,egnas, gkvbar,gkv,egkv, gl,el,ea,gabar , gst , gstbar , est
	GLOBAL hfinf, hftau,mftau ,mfinf ,hsinf, hstau,mstau ,msinf , ntau , ninf , mainf , hainf , matau , hatau , hstinf, hsttau,msttau ,mstinf

}
 
PARAMETER {
        gnafbar = 140 (mho/cm2)	<0,1e9>
	egnaf	= 58 (mV)
        gnasbar = 12 (mho/cm2)	<0,1e9>
	egnas	= 58 (mV)		
        gkvbar = 6 (mho/cm2)	<0,1e9>
	egkv	= -81 (mV)	
        gl = 0.40 (mho/cm2)	<0,1e9>
        el = -65 (mV)
        gstbar = 8.11 (mho/cm2)	<0,1e9>
        est = -81 (mV)		
        gabar = 58.1 (mho/cm2)	<0,1e9>
        ea = -81 (mV)	        	
}
 
STATE {
        mf msl hf hs n ma ha hst mst
}
 
ASSIGNED {
        v (mV)
	celsius (degC)

	gnaf (mho/cm2)
        inaf (mA/cm2)
	gnas (mho/cm2)
        inas (mA/cm2)
	gkv (mho/cm2)
        ikv (mA/cm2)
        il (mA/cm2)
        gst (mho/cm2)
        ist (mA/cm2)
        ga (mA/cm2)
        ia (mA/cm2)
        mfinf hfinf msinf hsinf ninf hainf mainf hstinf mstinf 
	hftau (ms) hstau (ms) ntau (ms) mftau (ms) mstau (ms) matau (ms) hatau (ms) msttau (ms) hsttau (ms)
}
 
LOCAL mfexp, hfexp, msexp, hsexp, mstexp, hstexp, nexp , haexp , maexp       
 
? currents
BREAKPOINT {
        SOLVE states METHOD cnexp
	
        gnaf = gnafbar*mf*mf*mf*hf
	inaf = gnaf*(v - egnaf)
        gnas = gnasbar*msl*msl*msl*hs
	inas = gnas*(v - egnas)
        gkv = gkvbar*(n^4)
        ikv = gkv*(v - egkv)
        ga = gabar*(ma^3)*ha
	ia  = ga*(v - ea)     
        gst = gstbar*(mst^3)*hst
        ist = gst*(v - est)
        il = gl*(v - el)
        
}
 
 
INITIAL {
	rates(v)
	mf = mfinf
	hf = hfinf
        msl = msinf
	hs = hsinf
	n = ninf
	ha = hainf
	ma = mainf
	hst = hstinf
	mst = mstinf        
}

? states
DERIVATIVE states {  
        rates(v)
        hf' = (hfinf-hf)/hftau
        mf' = (mfinf -mf)/mftau
        hs' = (hsinf-hs)/hstau
        msl' = (msinf -msl)/mstau
        n' = (ninf-n)/ntau
        ha' = (hainf-ha)/hatau
        ma' = (mainf -ma)/matau       
        hst' = (hstinf-hst)/hsttau
        mst' = (mstinf -mst)/msttau
}
 
LOCAL q10 , q11


? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                          :Call once from HOC to initialize inf at resting v.
		      
        LOCAL  alpha, beta, sum ,vd
        TABLE mfinf, hfinf, hftau, mftau, ninf, ntau, msinf, hsinf, hstinf, mstinf, hsttau, msttau   FROM -120 TO 40 WITH 641 :DEPEND celsius FROM -100 TO 100 WITH 200


UNITSOFF
        q10 = 3^((22.0 - 36)/10.0)
        q11 = 3^((36.0 - 23.5)/10.0)
        

               :"mf" sodium activation system
        vd = v 
        alpha = 0
        beta =  0
        sum = alpha + beta
        mfinf = 1.0 / (1 + exp((-30.1 - vd) / 6.65))
        mftau = (0.83 - 0.093) / (1 + exp((vd + 20.3) / 6.45)) + 0.093

                :"hf" sodium inactivation system
        vd = v 
        alpha = 0
        beta = 0
        sum = 0
	hftau = (1.66 - 0.12) / (1 + exp((vd + 8.03) / 8.69)) + 0.12
        hfinf = 1.0 / (1 + exp((vd + 51.4) / 5.9 ))

               :"msl" sodium activation system
        vd = v 
        alpha = 0
        beta =  0
        sum = alpha + beta
        msinf = 1.0 / (1 + exp((-30.1 - vd) / 6.65))
        mstau = (0.83 - 0.093) / (1 + exp((vd + 20.3) / 6.45)) + 0.093

                :"hs" sodium inactivation system
        vd = v 
        alpha = 0
        beta = 0
        sum = 0
	hstau = (12.24 - 1.9) / (1 + exp((vd + 32.6) / 8.0)) + 1.9
        hsinf = 1.0 / (1 + exp((vd + 51.4) / 5.9 ))

                :"n" Transient KV (delayed rectifier) current in honey bee from Wustenberg et al 2004
        vd = v 
        alpha = 0 
        beta = 0
	sum = 0
        ntau = (3.53 - 1.85) / (1 + exp((vd - 45) / 13.71)) + 1.85
        ninf = 1.0 / (1 + exp((-37.6 - vd) / 27.24))


                :"ha" Transient KA current in honey bee from Wustenberg et al 2004
        vd = v 
        alpha = 0
        beta = 0 
        sum = 0 
	hatau = (90 - 2.5) / ((1 + exp(- (vd + 60) / 25.0)) * (1 + exp((vd + 62) / 16.0))) + 2.5
        hainf = 1.0 / ( 1 + exp( ( vd + 74.7 ) / 7 ) )
        

                :"ma" Transient KA current in honey bee from Wustenberg et al 2004
        vd = v 
        alpha = 0
        beta = 0 
        sum = 0
	matau = (1.65 - 0.35) / ((1 + exp(- (vd + 70) / 4.0)) * (1 + exp((vd + 20) / 12.0))) + 0.35
        mainf = 1.0 / (1 + exp((-20.1 - vd)/16.1))


                :"hst" Slow transient KST current in honey bee from Wustenberg et al 2004
        vd = v 
        alpha = 0
        beta = 0 
        sum = 0 
	hsttau = (200.0 - 150.0) / (1 + exp((vd - 52) / 15.0)) + 150
        hstinf = 1.0 / ( 1 + exp( ( vd + 74.7 ) / 7 ) )
        

                :"mst" Slow transient KST current in honey bee from Wustenberg et al 2004
        vd = v 
        alpha = 0
        beta = 0 
        sum = 0
	msttau = (5.0 - 0.5) / (1 + exp((vd - 20) / 20.0)) + 0.5
        mstinf = 1.0 / (1 + exp((-20.1 - vd)/16.1))
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}

 
 
UNITSON
