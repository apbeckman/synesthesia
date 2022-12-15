//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 

const float det=.001;
const float maxdist=50.;
float l=0.;
mat2 rotm;
vec3 basecol=vec3(1.25,.75,1.);

mat3 lookat(vec3 dir, vec3 up){
    dir=normalize(dir);vec3 rt=normalize(cross(dir,normalize(up)));
    return mat3(rt,cross(rt,dir),dir);
}

mat2 rot2D(float a) {
    a=radians(a);
	float s=sin(a);
    float c=cos(a);
    return mat2(c,s,-s,c);
}


float de(vec3 p) {
    p=abs(5.-mod(p+5.,10.));
    float md=110.;
    float s=1.25;
    float sc=1.;
    vec3 pc;
	vec3 mp=vec3(95.);
    float rot=sin(smoothTimeC * .025)*17.5+ (cos(smoothTimeC * .025)*17.5);
    for (int i=0; i<10; i++) {
        p.xy=abs(p.xy); 
        p=p*s-1./sc;
        sc*=s;
        p.xz = _rotate(p.xz, RotateXY.x*PI);
        p.yz = _rotate(p.yz, RotateXY.y*PI);

        p.xz*=rotm+(RotateXY.x*0.125);
        p.yz*=rot2D((30.)+rot);
        float d=length(p.xz+sin(p.y)*.5)-.2/sc;
		mp=min(mp,abs(p));
        if (d<md) {
        	md=d;
			pc=p;
        }
    }
    l=mod(pc.y*.05-(smoothTimeC*0.75)*.1,.75)*2.;
    return md/sc;
}

vec3 march(vec3 from, vec3 dir) {
	vec3 p, col=vec3(0.);
    float totdist=0., d;
    for (int i=0; i<90; i++) {
    	p=from+totdist*dir;
        d=de(p);
    	totdist+=max(det,d);
        if (totdist>maxdist||length(col)>.3) break;
        col+=max(0.,det-d)*l;
    }
	col=.0+col*2.5*vec3(2.,3.,1.);
    return col;
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 uv=fragCoord/RENDERSIZE.xy-.5;
    uv.xy += FOV*_uvc*PI;
    uv.x*=RENDERSIZE.x/RENDERSIZE.y;
    uv += uv*_uvc*PI*PI*Warp;
    uv.x *=1.0+0.5*Warp.x;
    rotm=rot2D(-90.);
    vec3 dir=normalize(vec3(uv,.7));
    vec3 from=vec3(1.,2.,-5.);
        from.yz = _rotate(from.yz, Orbit.y*PI);

    from.xz = _rotate(from.xz, Orbit.x*PI);
    from.xz*=rot2D(smoothTime);
    from.yz*=rot2D(smoothTime);
    dir=lookat(-from,vec3(.5,1.,0.))*dir;
	vec3 col=march(from, dir);   
	col=mix(vec3(1.),col,min(1.,smoothTimeC*2.12));
  //  col=min(col,1.-smoothstep(.85,1.,abs(1.-mod(uv.y*1200.,2.)))*.4);
    fragColor = vec4(col,1.0);
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}