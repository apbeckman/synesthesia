

// created by florian berger (flockaroo) - 2023
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
//
//  lada taigaperfect loop - genuary1
//

#define LOOPING

//#define RandTex iChannel2
 
////////////////////////
//// quaternions, sdf's, helper funcs
////////////////////////

#define PI  3.14159265359
#define PI2 6.28318530718
#define PIH 1.57079632679

#define ROTM(ang) mat2(cos(ang-vec2(0,PIH)),-sin(ang-vec2(0,PIH)))

//#define ENABLE_MATERIALS 
#ifdef ENABLE_MATERIALS
#define SET_PREV_MAT(x) if(abs(d-d_mat)>.0001) mat=(x); d_mat=d;
#else
#define SET_PREV_MAT(x) 
#endif


vec3 rotZ(float ang,vec3 v) { return vec3(ROTM(ang)*v.xy,v.z); }

vec2 uvSmooth(vec2 uv,vec2 res)
{
    // no interpolation
    //return uv;
    // sinus interpolation
    return uv+1.*sin(uv*res*PI2)/(res*PI2);
    // iq's polynomial interpolation
    vec2 f = fract(uv*res);
    return (uv*res+.5-f+3.*f*f-2.0*f*f*f)/res;
}

vec4 inverseQuat(vec4 q)
{
    //return vec4(-q.xyz,q.w)/length(q);
    // if already normalized this is enough
    return vec4(-q.xyz,q.w);
}

vec4 multQuat(vec4 a, vec4 b)
{
    return vec4(cross(a.xyz,b.xyz) + a.xyz*b.w + b.xyz*a.w, a.w*b.w - dot(a.xyz,b.xyz));
}

vec3 transformVecByQuat( vec3 v, vec4 q )
{
    return (v + 2.0 * cross( q.xyz, cross( q.xyz, v ) + q.w*v ));
}

vec4 angVec2Quat(vec3 ang)
{
    float lang=length(ang);
    return vec4(ang/lang,1) * sin(vec2(lang*.5)+vec2(0,PI2*.25)).xxxy;
}

vec4 axAng2Quat(vec3 ax, float ang)
{
    return vec4(normalize(ax),1)*sin(vec2(ang*.5)+vec2(0,PI2*.25)).xxxy;
}

// iq's sdf primitives
float distBox( vec3 p, vec3 halfSize)
{
    vec3 q = abs(p) - halfSize;
    return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float distBoxR( vec3 p, vec3 halfSize, float r) { return distBox( p, halfSize-r ) - r ; }

float distCyl( vec3 p, float r, float h )
{
  vec2 d = vec2( length(p.xy)-r, abs(p.z) - h*.5 );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float distCylR( vec3 p, float r, float h, float R )
{
  vec2 d = vec2( length(p.xy)-(r-R), abs(p.z) - (h*.5-R) );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0))-R;
}

float distTorus(vec3 p, float R, float r)
{
    return length(p-vec3(normalize(p.xy),0)*R)-r;
}

float dDirLine(vec3 p, vec3 c, vec3 dir, float l)
{
    p-=c;
    dir=normalize(dir);
    float dp=dot(p,dir);
    //return length(p-dp*dir);
    return max(max(length(p-dp*dir),-dp),dp-l);
}

// iq's exponantial smooth-min func
float smin( float a, float b, float k )
{
    k=3./k;
    float res = exp2( -k*a ) + exp2( -k*b );
    return -log2( res )/k;
}

// iq's polynomial smooth-min func
float smin_( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

// flatness: // 0->sphere, 100-> nearly cylindric
float distTire(vec3 p, float r, float w, float h, float flatness)
{
    float l=length(p.xy);
    //p=abs(p);
    float d=1000.;
    // outer sphere
    float rfl=r*(1.+flatness);
    d=min(d,length(vec2(l+rfl-r,p.z))-rfl);
    float rz=-(rfl-r)+sqrt(rfl*rfl-p.z*p.z);
    //d=min(d,l-rz);
    float ang = atan(p.x,p.y);
    p.z+=cos(ang*64.)*w*.01*smoothstep(.87*r,1.*r,l);
    // main torus
    d=max(d,length(vec2(l-r+h*.5,p.z))-w*.5);
    //d=max(d,-l+r*.61);
    float w_l=sqrt(w*w-h*h); // w_laufflaeche
    float dz=.243*w_l;
    float zfr=mod(p.z,dz);
    float z=p.z-zfr+dz*.5;
    // rillen
    d=max(d,-length(vec2(l-rz,p.z-z))+dz*.2);
    // rim radius
    d=max(d,-(l-(r-h)));
    return d;
}

float distRim(vec3 p, float r, float w, float sh)  // outer rim radius, rim width;
{
    vec3 p0=p;
    p.z=abs(p.z);
    float ang0 = atan(p.y,p.x);

    float d=1000.,d2,d3;
    float dmain=length(p.xy)-r-sh;
    float dplane=-p.z+w*.5;
    d=-smin_(dplane,-dmain,.005);
    
    d2=-smin_(-(dmain+.005),-(dplane-.005),.01);
    d=-smin_(-d,d2,.01);
    
    d2=-smin_(-(dmain+.02),-(dplane-.045),.01);
    d=-smin_(-d,d2,.01);
    
    d2=dmain+.04;
    d=-smin_(-d,d2,.01);
    
    dplane=-p0.z+w*.5;
    d2=-smin_(-(dmain+.04),dplane+.015,.1);
    float c5=cos(ang0*5.);
    c5=-c5*.5+.5;
    c5*=c5;
    c5*=c5;
    c5=1.-c5;
    c5=mix(c5,0.,1.-clamp(((dmain+r+sh)-.045)/.03,0.,1.));
    d3=-smin_(-(dmain+.115),-(dplane-.01*(.8+.5*c5)),.04);
    d2=abs(-smin_(-d2,d3,.01))-.0015;
    d2=max(d2,-p0.z-.01);
    d=min(d,d2);
    
    float mang,ang;
    float dang;

    // rim holes
    dang=PI2/12.;
    mang=mod(ang0,dang);
    ang=ang0-mang+dang*.5;
    vec2 cs=cos(ang-vec2(0,PIH));
    d=max(d,-distCyl(p-vec3(r*.65*cs*(1.+3.*p.z),0.),.085*r,w*1.5));
    
    p=p0-vec3(0,0,.07);
    
    // screw holes
    dang=PI2/5.;
    mang=mod(ang0+dang*.5,dang);
    ang=ang0-mang+dang*.5;
    d2=distCyl(p-vec3(r*.3*cos(ang-vec2(0,PIH)),w*.05),.016,w*.19);
    // screws
    d=min(d, d2+.005);

    // axle
    d=min(d, distCyl(p,.045-.01,w*.25-.01)-.01);
    return d;
}

float distWheelDim(vec3 p, float w_mm, float h_perc, float rimD_inch, float shoulder_mm, float flatness)
{
    float w=w_mm*.001;
    float h=w*h_perc/100.;
    float d=10000.,d2;
    float rrim=rimD_inch*.5*.0254;
    d2=distTire(p, rrim+h, w, h, flatness );
    d=min(d,d2);
    float rimw=sqrt(w*w-h*h)+shoulder_mm*.001*2.5;
    d2=distRim(p, rrim, rimw, shoulder_mm*.001 );
    d=min(d,d2);
    return d;
}


bool intersectBox(vec3 p, vec3 dir, vec3 size)
{
    size*=.5*sign(dir);

    vec3 vmin = (-size-p)/dir;
    vec3 vmax = ( size-p)/dir;
    
    float tmin=vmin.x, tmax=vmax.x;
    
    if ((tmin > vmax.y) || (vmin.y > tmax)) return false; 
    tmin=max(tmin,vmin.y);
    tmax=min(tmax,vmax.y);
 
    if ((tmin > vmax.z) || (vmin.z > tmax)) return false; 
    tmin=max(tmin,vmin.z);
    tmax=min(tmax,vmax.z);
 
    return true; 
}
/*
vec4 getRand(vec2 coord)
{
    vec4 c=vec4(0);
    c+=texture(RandTex,coord+.003*TIME);
    c+=texture(RandTex,coord/2.+.003*TIME)*2.;
    c+=texture(RandTex,coord/4.+.003*TIME)*4.;
    c+=texture(RandTex,coord/8.+.003*TIME)*8.;
    return c/(1.+2.+4.+8.);
}
*/
#define FloorZ -.66
//#define HomePos vec3(0,0,-FloorZ*1.5)
//#define CamDist0 18.

// environment just a sky and some dark floor
vec4 myenv(vec3 pos, vec3 dir, float period_)
{
//#ifdef CUBEMAP    
//    return texture(iChannel1,dir.yzx);
//#endif
    vec3 sun = normalize(vec3(1,1,1));
    vec3 skyPos=pos+dir/abs(dir.z)*(120.-pos.z);
    //float cloudPat=(1.+.4*(getRand(skyPos.xy*.0002).x-.5));
    float cloudPat=(1.+.4);
    vec3 colHor=vec3(.3,.4,.5)+.4;
    float dirl=dot(dir,sun);
    vec3 clouds=mix(vec3(1.)*(1.-2.*dirl),vec3(.8,1.,1.2),cloudPat);
    vec3 colSky=mix(vec3(1.5,.75,0.)*3.,clouds,clamp(7.*dir.z,0.,1.));
    //colSky=mix(colSky,vec3(1),cloudPat);
    //colSky*=mix(1.,cloudPat,dir.z*5.);
    vec3 colFloor=vec3(.45);
    
    vec3 col=mix(colSky,colFloor,1.-smoothstep(-.01,.01,dir.z));
    col=mix(colHor,col,clamp(abs(dir.z*5.)-.1,0.,1.));
    
    col*=.9;
    
    //float sunang=acos(dot(dir,sun));
    float sunang=atan(length(cross(dir,sun)),dot(dir,sun));
    col+=15.*clamp(2.*exp(-sunang/.02),0.,1.);
    col+=2.*clamp(2.*exp(-sunang/.20),0.,1.);
    
    return vec4(col,1);
}



#define Res vec2(RENDERSIZE.xy)
//#define Res0 vec2(textureSize(iChannel0,0))
//#define Res1 vec2(textureSize(iChannel1,0))

#define enable_glass true

vec3 Delta=vec3(-.769,-1.073,-.669);

#define LoopNumFrames 240

float distCar(vec3 p)   
{
    p+=Delta;
    p=transformVecByQuat(p,axAng2Quat(vec3(1,0,0),.023));
    float d=100000.,d2;
    vec3 p00=p;
    p.x=abs(p.x);
    vec3 p0=p;
    vec3 psq=p*p;

    #ifdef ENABLE_MATERIALS 
    float d_mat=1001., mat=-1.;
    #endif
    SET_PREV_MAT(BG);
    
    vec3 frontWheelPos=vec3(1.44*.5,1.25,-.45);
    vec3 rearWheelPos=vec3(1.42*.5,-1.05,-.5);
    
    //d=min(d,distBox(p-vec3(0,0,-.05),vec3(1.65,3.45,.7)*.5));
    //d=min(d,distBox(p-vec3(0,-.5,.55),vec3(1.65,2.4,.55)*.5));
    bool front = p0.y>0.; 
    vec3 pwheel = front?frontWheelPos:rearWheelPos;
    // -- 15 ----- wheel cases (precalc) ------------
    float dWheelcases=distCylR((p-pwheel+vec3(0,0,-.1-(front?-.018:.0))).zyx,.4,.7,.1);
    
    float wcext=.01+(p.z+.5)*.1;

    // -- 0 ---------- main box -------------
    d=min(d,
          distBoxR(p-vec3(0,0,.2),
                   vec3(1.63-(psq.y*psq.y*.1+1.)*psq.y*.01-psq.z*.521-smoothstep(0.25,0.35,p.z)*(p.z-.25)*.0
                        +wcext*(1.-smoothstep(wcext*.3,wcext,dWheelcases))
//                        +wcext*exp(-dWheelcases*dWheelcases/wcext/wcext)
                        ,
                        3.6- min(psq.z*psq.z*20.*step(0.,p.y),.1) - (psq.z+psq.x)*.1  -step(0.,-p.y)*smoothstep(0.,.2,p.z)*max(p.z-.2,0.)*1.1,
                        1.3-psq.y*psq.y*psq.y*.015*(.5+.5*step(0.,-p.z))-step(0.,p.z)*(psq.x*.2+psq.y*.035+p.y*.04))*.5,
                        .03-min(p.y*.03,0.)));
    float dmainBox=d;
    
    // -- 1 ---------- hood/windshield cut -------------
    float dwin =dot(p-vec3(0,1.5-psq.x*(.3-(p.z-.8)*.4),0.7)*.5,normalize(vec3(0,1.,1)));
    float dhood=dot(p-vec3(0,1.5-psq.x*(2.-(p.y-.75)*1.),0.74-.1*(p.y-.75)*(p.y-.75))*.5,normalize(vec3(0,.07,1)));
    //dwin=10000.;
    d=-smin_(-d,-smin_(dwin,dhood,.03),max(.001,.05-.0*(p.y*.01)));
    //d=-min(-d,-min(dwin,dhood));

    float dwin2 =dot(p-vec3(0,1.64-psq.x*.3,0.7)*.5,normalize(vec3(0,1.5,1)));
    float dhood2=dot(p-vec3(0,1.5,0.36)*.5,normalize(vec3(0,.04,1)));
    d2=-min(-d,-min(dwin2,dhood2));
    d+=exp(-abs(d2)/.006)*.006;
    
    p=p-vec3(0,0,-.01*psq.y);

    // -------- absatz hood ---------------    
    d2=distBoxR(p-vec3(0.,1.37,.3),vec3(.9,1,.5)*.5,.07);
    d+=.015*(1.-smoothstep(-.015,.015,d2))*clamp(1.7-p0.y,0.,1.);
    
    // -- 2 ---------- side versenkung windows --------------
    d2=distBoxR(p-vec3(.9-p.z*.2,0,.57),vec3(.1,3.25-p.z*.75,.59)*.5+.02,.15)-.02;
    d=-smin_(-d,d2,.01);

    SET_PREV_MAT(CARBODY);
    
    float dwincut=1000.;
    // -- 3 ---------- side win rear -----------
    p=p-vec3(0,-.88,.55);
    d2=distBoxR(p,vec3(3.,.82-.8*(p.z)*step(0.,-p.y),.37)*.5,.063-p.y*.05);
    //d+=clamp(-d2,0.,.01);
    //d+=smoothstep(0.,1.,-d2/.01)*.01;
    ////d+=exp(-abs(d2)/.006)*.006;
#ifdef ENABLE_MATERIALS 
    //if(enable_glass) { SET_PREV_MAT(GLASS); }
#endif
    dwincut=min(dwincut,d2);
    p=p0;

    // -- 4 ---------- side stripe --------------
    d2=distBox(p-vec3(0,0,.16-.01*(p.y+.4)*(p.y+.4)),vec3(3.,3.3,.07)*.5);
    //d+=clamp(-d2,0.,.01);
    d+=smoothstep(0.,1.,-d2/.01)*.01;

    // -- 5 ---------- door -------------
    d2=distBoxR(p-vec3(0,.1,.235),vec3(3.,1.03,1.09)*.5,.07-.1*(p.z+.3));
    d2=-smin_(-(dwin+.03),-d2,.07);
    d+=exp(-abs(d2)/.006)*.01;
    SET_PREV_MAT(CARBODY);

    // -- 6 ---------- side win front -----------
    d2=-smin_(-d2-.05,p.z-.36,.03);
    //d+=clamp(-d2,0.,.01);
    //d+=smoothstep(0.,1.,-d2/.01)*.01;
    ////d+=exp(-abs(d2)/.006)*.006;
    dwincut=min(dwincut,d2);
    
    // -- 7 ---------- front window -----------
    p=p0-vec3(0,0,.53);
    //d2=distBoxR(p,vec3(1.4,5.,.37)*.5,.05);
    d2=dmainBox+.065-psq.x*.1*p.z;
    d2=-smin(-d2,dhood-.015,.05);
    //d+=clamp(-d2,0.,.01);
    //d+=smoothstep(0.,1.,-d2/.015)*.01;
    dwincut=min(dwincut,d2);
#ifdef ENABLE_MATERIALS 
    //if(enable_glass) { SET_PREV_MAT(GLASS); }
#endif

    // -- 8 ---------- rear window -----------
    p=p0-vec3(0,-1.7,.5);
    d2=distBoxR(p,vec3(1.15-p.z*.3,1.,.35)*.48,.07);
    d+=smoothstep(0.,1.,-(d2-.03+.03*p.z)/.01)*.01;
    SET_PREV_MAT(CARBODY);
    //d+=smoothstep(0.,1.,-d2/.01)*.01;
    dwincut=min(dwincut,d2);
    //d=max(d,-d2);
#ifdef ENABLE_MATERIALS 
    //if(enable_glass) { SET_PREV_MAT(GLASS); }
#endif

    // ------- bottom absatz ------------
    d2=p0.z+.345;
    d+=.02*min(exp2(-d2/.02),1.);

    // -- 15 ----- wheel cases (apply) ------------
    d=-smin_(-d,dWheelcases,.005);
    
    SET_PREV_MAT(p0.z<-.345?GUMMI:CARBODY);

    // ------------- cut out interior, cutout windows, add windows + window lips
    float d_inner=d+.06;
    d=max(d,-d_inner); SET_PREV_MAT(BLACKPLASTIC); 
    d=-smin_(-d,dwincut,.01);
    SET_PREV_MAT(CARBODY); 
    d=min(d,length(vec2(dwincut,d_inner-.045))-.01);
    SET_PREV_MAT(GUMMI); 
    d=min(d,length(vec2(dwincut,d_inner-.045)+vec2(1,-1)*.008)-.0025);
    SET_PREV_MAT(CHROME); 
    if(enable_glass) { d=min(d,d_inner-.045); SET_PREV_MAT(GLASS); }
    
    // -- 9 ---------- rear door -----------
    p=p0-vec3(0,-1.7,.25);
    d2=distBoxR(p,vec3(1.35-psq.z*.35,.5+p0.z*.5,.97)*.5,.07);
    d+=exp(-abs(d2)/.006)*.01;

    // -- 11 ---------- rear license plate box -----------
    p=p0-vec3(0,-1.885,.04);
    d2=distBoxR(p,vec3(.65+p.z*.3,.2,.2)*.5,.03);
    d=max(d,-d2);
    SET_PREV_MAT(CARBODY);
    p=p0-vec3(0,-1.775,.14);
    d2=distBoxR(p,vec3(.68,.05,.04)*.5,.005);
    d=min(d,d2);

    // -- 10 ---------- rear blinker -----------
    p=p0-vec3(.68,-1.8,.0);
    d2=distBoxR(p,vec3(.2-p.z*.2*step(0.,-p.x),.3,.3)*.5,.03);
    d=min(d,d2);
    d=max(d,dmainBox-.01);
    SET_PREV_MAT(BLACKPLASTIC);
    d-=smoothstep(.007,.01,-d2)*.001;
    SET_PREV_MAT(p.z<0.06?(p.x>0.0?(p.z>-.05?GLASS:ORANGEGLASS):REDGLASS):REDGLASS);
    
    // -- 12 ---------- front blinker -----------
    p=p0-vec3(.6,1.75,.15);
    d2=distBoxR(p,vec3(.25,.095,.12)*.5,.03);
    d=max(d,-d2);
    SET_PREV_MAT(BLACKPLASTIC);
    d=min(d,-smin_(-d2-.025,-p.y+.01,.01));
    SET_PREV_MAT(p.x<0.03?ORANGEGLASS:GLASS);

    // -- 13 ---------- grill -----------
    p=p0-vec3(.0,1.78-psq.x*.05,-.05);
    d2=distBoxR(p,vec3(1.48+.02-p.y*.3,.2,.265-psq.x*.03+.02-p.y*.3)*.5,.06-p.z*.15);
    d=max(d,-d2);
    // ----- einbuchtung rund um scheinwerfer -----
    float d2b=-smin_(-d2,p.x-.48,.02);
    vec2 pl=(p-vec3(.61,0,0)).xz; float lpl=length(pl);
    float yo=mix(.02,max(-d2b-.016,0.)*.5,1.-exp2((-lpl+.1)/.0075));
    d2=-smin_(-d2,-.1+lpl,.005);
    d2=-smin_(-d2,-abs(p.y-.01)+.035-yo,.005);
    // ----- cooling slits -----
    vec3 pi=vec3(0,0,(clamp(floor(p.z/.029),-4.,3.)+.5)*.029);
    d2=-smin_(-d2,distBox(p-pi,vec3(.96,.2,.02)*.5),.01);
    // ----- lada logo -------
    p.y-=.03;
    d2=min(d2,distBoxR(p,vec3(.04+p.z*.1,.02,.07),.005));
    d=min(d,d2);
    SET_PREV_MAT(BLACKPLASTIC);

    // ----- scheinwerfer -----
    float lsph=length(p-vec3(.61,-.273-.03,0));
    d=min(d,max(max(lsph-.3,-p0.y+1.65),lpl-.09));
    SET_PREV_MAT(CHROME);

    // -- 14 ---------- bumpers -----------
    p=p0; p.y=abs(p.y);
    p=p-vec3(0,1.87-.03*psq.x,-.25);
    d2=distBoxR(p,vec3(1.62,.1,.1)*.5,.02);
    d=min(d,d2);
    SET_PREV_MAT(CHROME);
    p-=vec3(.76-p.y*.15,-.065,0);
    d2=distBoxR(p,vec3(.113,.25,.113)*.5,.02);
    d2=max(d2,dot(p,vec3(-1,-.4,0))-.03);
    d=min(d,d2);
    SET_PREV_MAT(GUMMI);
    
    // ------- wheels, axes ------------
    //d=max(d,-distCylR((p-rearWheelPos+vec3(0,0,-.1)).zyx,.4,.7,.1));
    // ------- axes ------
    #if 0
    p=(p0-pwheel*vec3(0,1,1));
#ifdef USE_SIMDATA
    float leftSgn=sign(p00.x);
    float rear=front?0.:1.;
    vec3 wo=vec3(mix(WheelDistF,WheelDistR,rear)*.5,0,0);

    // wheel offsets
    vec4 qf=axAng2Quat(vec3(0,1,0),(WheelOffsFL-WheelOffsFR)/WheelDistF*leftSgn);
    vec4 qr=axAng2Quat(vec3(0,1,0),(WheelOffsRL-WheelOffsRR)/WheelDistR*leftSgn);
    vec4 axQuat=front?qf:qr;
    float axOffs=mix(WheelOffsFR+WheelOffsFL,WheelOffsRR+WheelOffsRL,rear)*.5;

    p=transformVecByQuat(p,axQuat);
    p.z-=axOffs;
#endif
    d=min(d,distCylR(p.zyx,.12*(1.-.6*smoothstep(0.05,.2,p.x)),1.44,.05));
    #endif
    p=(p0-pwheel);

    vec3 wo=vec3(0,0,0);
#ifdef USE_SIMDATA
    float leftSgn=sign(p00.x);
    float rear=front?0.:1.;
    wo=vec3(mix(WheelDistF,WheelDistR,rear)*.5,0,0);

    // wheel offsets
    vec4 qf=axAng2Quat(vec3(0,1,0),(WheelOffsFL-WheelOffsFR)/WheelDistF*leftSgn);
    vec4 qr=axAng2Quat(vec3(0,1,0),(WheelOffsRL-WheelOffsRR)/WheelDistR*leftSgn);
    vec4 axQuat=front?qf:qr;
    float axOffs=mix(WheelOffsFR+WheelOffsFL,WheelOffsRR+WheelOffsRL,rear)*.5;
    // wheel axis rot + offset
    p=transformVecByQuat(p+wo,axQuat)-wo;
    p.z-=axOffs;
#endif
    d=min(d,distCyl(p.zyx-vec3(0,0,-.72),.12*(1.-.6*smoothstep(0.05,.2,p.x+.72)),1.55));
    SET_PREV_MAT(CHASSIS);
#ifdef USE_SIMDATA
    // steering rotation of front wheels
    vec4 q=axAng2Quat(vec3(0,0,1),leftSgn*(1.-.1*leftSgn*sign(SteerAng))*-SteerAng*(front?1.:0.));
    p=transformVecByQuat(p+vec3(.1,0,0),q)-vec3(.1,0,0);
    // wheel rotations
    float rot=-WheelRot.x*.7;
    p=transformVecByQuat(p,axAng2Quat(vec3(1,0,0),rot));
#endif
#define PROPER_WHEELS
#ifndef PROPER_WHEELS
    d=min(d,distCylR(p.zyx,.35,.2,.05));
#else
    // newer viva ...broader tires
    p.yz*=ROTM(-float(FRAMECOUNT)/float(LoopNumFrames)*PI2*3.); 
    d=min(d, distWheelDim(p.yzx,185.,75.,16.,12.,.2));
    // old niva
    //d=min(d, distWheelDim(p.yzx,175.,80.,16.,12.,.2));
#endif
    SET_PREV_MAT(length(p.yz)<.215?CHROME:GUMMI);
    
#ifdef ENABLE_MATERIALS 
    return vec2(d,mat);
#else
    return d;
#endif
}


/*float dist(vec3 pos)
{
    float fact=pow(100.,fract(-float(FRAMECOUNT)/float(LoopNumFrames)));
    pos*=fact;
    float R=2.,r=.01;
    float d=100000.;
    
    d=min(d,distCar(pos));
    d=min(d,distCar((pos+Delta*00.)/100.)*100.);
    d=min(d,distCar((pos+Delta*0000.)/10000.)*10000.);
    
    return d/fact;
}*/

// iq's superb compile optimization - thanks a lot!!
float dist(vec3 pos)
{
    float fact=pow(100.,fract(-float(FRAMECOUNT)/float(LoopNumFrames)));
    pos*=fact;
    float R=2.,r=.01;
    float d=100000.;
    
    for( int i=min(0,int(FRAMECOUNT)); i<3; i++ )
    {
        float sca = pow(100.0,float(i));
        d=min(d,distCar((pos+Delta*00.)/sca)*sca);
    }
    
    return d/fact;
}

/*vec3 getGrad(vec3 p, float eps) 
{ 
    vec2 d=vec2(eps,0); 
    float d0=dist(p);
    return vec3(dist(p+d.xyy)-d0,dist(p+d.yxy)-d0,dist(p+d.yyx)-d0)/eps; 
}*/

// iq's superb compile optimization - thanks a lot!!
vec3 getGrad(vec3 p, float eps) 
{ 
    vec3 n = vec3(0.0);
    for( int i=min(0,int(FRAMECOUNT)); i<4; i++ )
    {
        vec3 e = 0.5773*(2.0*vec3((((i+3)>>1)&1),((i>>1)&1),(i&1))-1.0);
        n += e*dist(p+eps*e);
    }
    return n;
}

float march(inout vec3 pos, vec3 dir)
{
    float eps=.003;
    for(int i=0;i<120;i++)
    {
        float d=dist(pos);
        pos+=dir*d*.6;
        if (d<eps) return 1.;
    }
    return 0.;
}

#define MOUSE_PHI   (_mouse.x/Res.x*12.)
#define MOUSE_THETA (_mouse.y/Res.y*12.)

void getTrafo(inout vec3 pos, inout vec3 dir, vec2 fc)
{
    vec2 sc=(fc-Res*.5)/Res.x*2.;
    float tanFOVh=1./2.5;
    dir=normalize(vec3(0,0,-1./tanFOVh)+vec3(sc,0));
    pos=vec3(0,0,6.1);
    #ifdef LOOPING
    pos.xy*=ROTM(-1.057-float(int(FRAMECOUNT))/float(LoopNumFrames)*PI2);
    dir.xy*=ROTM(-1.057-float(int(FRAMECOUNT))/float(LoopNumFrames)*PI2);
    #endif
    float ph = MOUSE_PHI;
    float th = MOUSE_THETA;
    if (_mouse.x<1.) { 
        ph=3.6;
        #ifdef LOOPING
        ph=4.;
        #endif
        th=1.3;
    }
    pos.yz=ROTM(th)*pos.yz;
    dir.yz=ROTM(th)*dir.yz;
    pos.xy=ROTM(ph)*pos.xy;
    dir.xy=ROTM(ph)*dir.xy;
}

vec4 ovlCol( vec2 fragCoord )
{
    vec4 col=vec4(0);
    vec2 sc=(fragCoord-RENDERSIZE.xy*.5)/RENDERSIZE.x;
    float m2sc=.205;
    
    float ph04=fract(MOUSE_PHI/PI*2.*.25)*4.;
    float backFact =max(clamp(-3.+ph04,0.,1.),clamp(1.-ph04,0.,1.));
    float sideFact =min(clamp( 0.+ph04,0.,1.),clamp(2.-ph04,0.,1.));
    float frontFact=min(clamp(-1.+ph04,0.,1.),clamp(3.-ph04,0.,1.));
    float sideFact2=min(clamp(-2.+ph04,0.,1.),clamp(4.-ph04,0.,1.));
    sideFact=max(sideFact,sideFact2);
    float scFact=sideFact2>0.?1.:-1.;
    //col+=sideFact *texture(iChannel0,(ROTM(.023*scFact)*sc*vec2(scFact,1))/Res0*556.+.5+vec2(-.033,.1));
    //col+=frontFact*texture(iChannel1,sc/Res1*556.+.5+vec2(-.252,.1));
    //col+=backFact *texture(iChannel1,sc/Res1*556.+.5+vec2(.26,.1));
    col.xyz-=exp(-sc.x*sc.x/.001/.001);
    col.xyz-=exp(-sc.y*sc.y/.001/.001);
    vec2 sinsc=sin(sc/m2sc*2.*PI);
    col.yz-=.3*exp(-sinsc.x*sinsc.x/.03/.03);
    col.yz-=.3*exp(-sinsc.y*sinsc.y/.03/.03);
    col.xz-=exp(-(sc.x-3.74*.5*m2sc)*(sc.x-3.74*.5*m2sc)/.001/.001);
    col.xz-=exp(-(sc.x+3.74*.5*m2sc)*(sc.x+3.74*.5*m2sc)/.001/.001);
    col.xz-=exp(-(sc.y-1.64*.5*m2sc)*(sc.y-1.64*.5*m2sc)/.001/.001);
    col.xz-=exp(-(sc.y+1.64*.5*m2sc)*(sc.y+1.64*.5*m2sc)/.001/.001);
    return col;
}

uniform float ovlFade;

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec4 ovl=ovlCol( fragCoord );
    vec3 pos,dir;
    getTrafo(pos,dir,fragCoord);
    float m=march(pos,dir); //m=1.;
    vec3 left=normalize(cross(vec3(0,0,1),dir));
    vec3 up  =normalize(cross(left,dir));
    vec3 g=getGrad(pos,.0005);

    float lg=length(g);
    vec3 n=g/(lg+.00001);
    
    float ao=1.;
    float sc=1.7,scmin=.085;
    float df=dist(pos+n*scmin*.2)/(scmin*.2);
    for(int i=0;i<20;i++)
    {
        float ao2=clamp(dist(pos+n*sc)/(sc/df),0.,1.);
        ao*=mix(ao2,1.,
        1.-sqrt(sc)*.2
        );
        sc*=.7;
        if(sc<scmin) break;
    }
    ao*=m*.5+.5;
    
    fragColor.xyz=(n*.15+.85);
    
    vec3 R=reflect(dir,n);
    vec4 refl=myenv(vec3(0),R,1.);
    float fres=1.-abs(dot(R,n));
    fres*=fres*fres;
    fres=.1+.7*fres;
    fragColor=mix(fragColor,refl*(ao*.35+.65),fres);
    
    fragColor.xyz*=ao;
    fragColor=mix(clamp(fragColor,0.,1.),ovl,ovlFade);
    
    fragColor.w=1.;
	return fragColor; 
 } 




vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}