//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


#define R(p,a,t) mix(a*dot(p,a),p,cos(t))+sin(t)*cross(p,a)
#define H(h) (cos((h)*6.3+vec3(0,23,21))*.5+.5)

vec4 renderMainImage() {
	vec4 O = vec4(0.0);
	vec2 C = _xy;
    vec2 r= RENDERSIZE.xy;
    vec3 p,c=vec3(0),
    d=normalize(vec3((C-.5*r.xy)/r.y,.7));
    float i=0.,s,e,g=0.,t=TIME;
	for(;i++<120.;){
        p=R(g*d,normalize(H(smoothTimeC*.05)),g*.1);
        //rotating the camera using what I think is fabrice's ray rotating matrix - it's used in a lot of Shane's stuff, and I've 
        //copied+pasted it from Maze Lattice and modified it to fit a shader more times than I can count
        p.xy = _rotate(p.xy, lookXY.x*PI);//spins AB
        p.yz = _rotate(p.yz, lookXY.y*PI);//up/down AB
        p.z+=smoothTime*.5; //reactive camera movement - there are three of these uniforms in the script.js, smoothTime is bass, ' 'B is highs, and ' 'C is mids
        p=asin(.7*sin(p));
        s=2.5+sin(.5*smoothTime+3.*sin(t*1.)*0.125)*.5;
        for(int i=0;i++<(6+itera);p=p*e-vec3(3,2.5,3.5))
            p=abs(p+test*PI),
            p=p.x<p.y?p.zxy:p.zyx,
            s*=e=2.+iter;
        g+=e=abs(length(p.xz)-.3)/s+2e-5;
	    c+=mix(vec3(1),H(p.z*.5+smoothTime*.1),.4)*.02/exp(.5*i*i*e)*( 1.0+highhits);
	}
	c*=c*( 1.0+highhits);
    O=vec4(c,1);
   // O*=( 1.0+highhits);
	return O; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}