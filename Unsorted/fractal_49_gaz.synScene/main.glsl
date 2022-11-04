

#define R(p,a,r)mix(a*dot(p,a),p,cos(r))+sin(r)*cross(p,a)
#define H(h)(cos((h)*6.3+vec3(0,23,21))*.5+.5)
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 C = _xy;

    O=vec4(0);
    vec3 p,r=vec3(RENDERSIZE, 0),
    d=normalize(vec3((C-.5*r.xy)/r.y,2));  
    for(float i=0.,s,e,g=1.5;
        ++i<90.;
        O.xyz+=.1*mix(vec3(1),H(log(s)*.3),.8)*exp(-12.*i*i*e)
    )
    {
        p=g*d-vec3(-.2,.3,2.5);
        p=R(p,normalize(vec3(1,2.*sin(TIME*.1),3)),TIME*.2);
        s=5.;
        p=p/dot(p,p)+1.;
        for(int i=0;i++<8;p*=e)
            p=1.-abs(p-1.),
            s*=e=1.6/min(dot(p,p),1.5);
        g+=e=length(cross(p,vec3(.577)))/s-5e-4;
    }
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}