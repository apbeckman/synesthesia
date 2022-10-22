

#define R(p,a,r)mix(a*dot(p,a),p,cos(r))+sin(r)*cross(p,a)
#define H(h)(cos((h)*6.3+vec3(0,23,21))*.5+.5)
vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 C = _xy;

    O=vec4(0);
    vec3 p,r=vec3(RENDERSIZE,0),n=vec3(-.5,-.809,.309),
    d=normalize(vec3((C-.5*r.xy)/r.y,1));  
    for(float i=0.,e,g=0.;
        ++i<99.;
        O.xyz+=mix(vec3(1),H(dot(p,p)*.5),.7)*.05*exp(-.05*i*i*e)
    )
    {
        p=g*d;
        p.z-=11.;
        p=R(p,normalize(vec3(1,2,2)),smoothTime*.25);
        for(int j=0;j<7;j++)
            p.xy=abs(p.xy),
            p-=2.*min(0.,dot(p,n))*n;
        p.z=fract(log(p.z)-smoothTime*.25)-.5;
        g+=e=length(p.yz)-.01;
        // Dodecahedron Frame
        //g+=e=length(p.xz)-.01;
    }
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}