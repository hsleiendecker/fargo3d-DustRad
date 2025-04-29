if (PhysicalTime>0.0 && id2==1 && !(CPU_Rank==0 && j<NGHY)){

    St1 = 1/InvSt1;
    St2 = 1/InvSt2;

    //relative drift velocity 
    v01 = pow(cs[l] * cs[l]/(2*r*omega)*gamma*(St1-St2),2);
    //dvy = pow(velocities_input[1][l]-velocities_input[2][l],2);
    //if(dvy>v01) v01=dvy;
    //v11 = pow(cs[l] * cs[l]/(2*r*omega)*gamma*(St2/2),2);
    //this ^ alpha_2pop = |alpha-2*alpha| bc it's (St1-St2) & St2=St1/2

    //relative settling v: <vsett> = rho_s*a*omega^2*<z>/(rho_g_@<z> * v_th)
    //2.507 is sqrt(2*pi) and 2.91 = sqrt(2*pi)*e^(-1/pi)*sqrt(8/pi) 
    //v01 += pow(rho_s*(a1-a0)*omega*omega*2.507*H_1/(rho[0][l]/(2.91*H_1*cs[l])),2); 
    //v11 = rho_s*(a1-a0)*omega*omega*2.507*H_1/(rho[0][l]/(2.91*H_1*cs[l]));
    
    if(St1>0.5) St1=0.5;
    if(St2>0.5) St2=0.5;
    v01 += pow((H_0*St1-H_1*St2)*omega,2);
    //v11 += pow((H_1*St2/2)*omega,2);

    //brownian motion: eqn 1 Ormel+ 2007
    v01 += pow(cs[l]*sqrt(8/3.1415*mu*m_p*(a1*a1*a1+a0*a0*a0)/(4.189*rho_s*a1*a1*a1*a0*a0*a0)),2);
    //v11 += pow(cs[l]*sqrt(8/3.1415*mu*m_p*2/(4.189*rho_s*a1*a1*a1)),2);

    //turbulent gas motion from Ormel & Cuzzi 2007
    //repeat variables
    OmKinv = 1.0 / omega;
    Re = 0.5 * ALPHA * rho[0][l] * sig_h2 / (mu*m_p);
    ReInvSqrt = sqrt(1.0 / Re);
    vn = sqrt(ALPHA) * cs[l];
    vs = pow(Re, -0.25) * vn;
    ts = OmKinv * ReInvSqrt;
    vg2 = 1.5 * vn * vn;
    //solve for Large->Small and Large->Medium relative velocities
    StL = 1/InvSt2;
    StS = 1/InvSt1;
    StM = StL/2;
    epsLS = StS / StL;
    epsLM = StM / StL;
    tauL = StL * OmKinv;
    tauM = StM * OmKinv;
    tauS = StS * OmKinv;
    ys = c0 + c1*StL + c2*StL*StL + c3*StL*StL*StL;
    h1LS = (StL - StS) / (StL + StS) *
         (StL * yap1inv - StS*StS / (StS + ya * StL));
    h1LM = (StL - StM) / (StL + StM) *
         (StL * yap1inv - StM*StM / (StM + ya * StL));
    h2LS = 2.0 * (ya * StL - ReInvSqrt) + StL * yap1inv -
         StL*StL / (StL + ReInvSqrt) + StS*StS / (ya * StL + StS) -
         StS*StS / (StS + ReInvSqrt);
    h2LM = 2.0 * (ya * StL - ReInvSqrt) + StL * yap1inv -
         StL*StL / (StL + ReInvSqrt) + StM*StM / (ya * StL + StM) -
         StM*StM / (StM + ReInvSqrt);
    //determine which Stokes regime for which to calculate these drift velocities
    int which_turb = 0;
    if (tauL < 0.2 * ts) {
      v01 += 1.5 * pow((vs / ts * (tauL - tauS)), 2);
      //v11 += 1.5 * pow((vs / ts * (tauL - tauM)), 2);
    } 
    else if (tauL * ya < ts) {
      v01 += vg2 * (StL - StS) / (StL + StS) *
              (StL*StL / (StL + ReInvSqrt) - StS*StS / (StS + ReInvSqrt));
      //v11 += vg2 * (StL - StM) / (StL + StM) *
      //        (StL*StL / (StL + ReInvSqrt) - StM*StM / (StM + ReInvSqrt));
    } 
    else if (tauL < 5 * ts) {
      //vRelLS = sqrt(vg2 * (h1LS + h2LS));
      //vRelLM = sqrt(vg2 * (h1LM + h2LM));
      v01 += vg2 * (h1LS + h2LS);
      //v11 += vg2 * (h1LM + h2LM);
    } 
    else if (tauL < 0.2 * OmKinv) {
      v01 += vg2 * StL * (2 * ya - 1 - epsLS + 2 / (1 + epsLS) *
                      (yap1inv + pow(epsLS, 3) / (ya + epsLS)));
      //v11 += vg2 * StL * (2 * ya - 1 - epsLM + 2 / (1 + epsLM) *
      //                (yap1inv + pow(epsLM, 3) / (ya + epsLM)));
    } 
    else if (tauL < OmKinv) {
      v01 += vg2 * StL * (2 * ys - 1 - epsLS + 2 / (1 + epsLS) *
                      (1 / (1 + ys) + pow(epsLS, 3) / (ys + epsLS)));
      //v11 += vg2 * StL * (2 * ys - 1 - epsLM + 2 / (1 + epsLM) *
      //                (1 / (1 + ys) + pow(epsLM, 3) / (ys + epsLM)));
    } 
    else {
      v01 += vg2 * (2.0 + StL + StS) /(1.0 + StL + StS + StL * StS);
      //v11 += vg2 * (2.0 + StL + StM) /(1.0 + StL + StM + StL * StM);
      which_turb = 6;
    }

    //take sqrt to complete adding in quadrature
    v01 = sqrt(v01);
    //v11 = sqrt(v11);

    sig01 = 3.1415*(a1*a1+a0*a0);
    mass_sweep = rho[1][l]*rho[2][l]*sig01*v01
                 /(4.189*rho_s*a1*a1*a1*sqrt(6.283*(H_0*H_0+H_1*H_1))) * dt;
    //if(mass_sweep<(pow(exp(tau_grow_ph),3)-1)*rho[1][l]) mass_sweep = (pow(exp(tau_grow_ph),3)-1)*rho[1][l];
    //if(rho[2][l]/rho[1][l]*plaw>1)
    //sig11 = 6.283*a1*a1;
    //F = sqrt(2*H_1*H_1/(H_0*H_0+H_1*H_1))*sig01/sig11*v01/v11*pow(a1/a0,0.5*-(4.0-plaw));
    //F = sqrt(2*H_1*H_1/(H_0*H_0+H_1*H_1))*sig01/sig11*v01/v11*plaw;
    //mass_frag = F*rho[2][l]*rho[2][l]*sig11*v11
      //          /(4.189*rho_s*a1*a1*a1*sqrt(12.566*H_1*H_1)) * dt;
    //after simplifying, mass frag is simply:
    mass_frag = mass_sweep * rho[2][l]/rho[1][l]*plaw;

    rho[1][l] -= mass_sweep;
    rho[2][l] += mass_sweep;

    rho[2][l] -= mass_frag;
    rho[1][l] += mass_frag;

}
