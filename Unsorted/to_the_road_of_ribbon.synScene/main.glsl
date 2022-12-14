//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


float tunnel(vec3 p)
{
	return cos(p.x)+cos(p.y*1.5)+cos(p.z)+cos(p.y*20.)*.05;
}

float ribbon(vec3 p)
{

	return length(max(abs(p-vec3(cos(p.z*1.5)*.3,-.5+cos(p.z)*.2,.0))-vec3(.125,.02,smoothTime+3.),vec3(.0)));
}

float scene(vec3 p)
{
	return min(tunnel(p),ribbon(p));
}

vec3 getNormal(vec3 p)
{
	vec3 eps=vec3(.1,0,0);
	return normalize(vec3(scene(p+eps.xyy),scene(p+eps.yxy),scene(p+eps.yyx)));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	vec2 v = -1.0 + 2.0 * fragCoord.xy / RENDERSIZE.xy;
	v += (_uvc*PI*FOV);
	v.x *= RENDERSIZE.x/RENDERSIZE.y;
 
	vec4 color = vec4(0.0);
	vec3 org   = vec3(sin(TIME*0.2)*.25,cos(TIME*.125)*.125+.25,smoothTime);
	vec3 dir   = normalize(vec3(v.x*1.6,v.y,1.0));
	dir.xy *= 1.0+((_uvc+dir.xy)-1.)*Warp;
	dir.xy = _rotate(dir.xy, Rotate*PI);

	dir.xz = _rotate(dir.xz, LookXY.x*PI);
	dir.yz = _rotate(dir.yz, LookXY.t*PI);
	
	vec3 p     = org,pp;
	float d    = .0;

	//First raymarching
	for(int i=0;i<64;i++)
	{
	  	d = scene(p);
		p += d*dir;
	}
	pp = p;
	float f=length(p-org)*0.02;

	//Second raymarching (reflection)
	dir=reflect(dir,getNormal(p));
	p+=dir;
	for(int i=0;i<48;i++)
	{
		d = scene(p);
	 	p += d*dir;
	}
	color = max(dot(getNormal(p),vec3(.1,.1,.0)), .0) + vec4(.3,cos(smoothTimeB*.125)*.5+.5,sin(smoothTimeB*.125)*.5+.5,1.)*min(length(p-org)*.04, 1.);

	//Ribbon Color
	if(tunnel(pp)>ribbon(pp))
		color = mix(color, vec4(cos(smoothTimeB*.125)*.5+.5,cos(smoothTimeB*.125)*.5+.5,sin(smoothTimeB*.125)*.5+.5,1.),.3);

	//Final Color
	vec4 fcolor = ((color+vec4(f))+(1.-min(pp.y+1.9,1.))*vec4(1.,.8,.7,1.))*min(TIME*.5,1.);
	fragColor = vec4(fcolor.xyz,1.0);
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}