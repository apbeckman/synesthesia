vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
// Created by S. Guillitte 2015

float zoom=1.;

vec2 cmul( vec2 a, vec2 b )  { return vec2( a.x*b.x - a.y*b.y, a.x*b.y + a.y*b.x ); }
vec2 csqr( vec2 a )  { return vec2( a.x*a.x - a.y*a.y, 2.*a.x*a.y  ); }


mat2 rot(float a) {
	return mat2(cos(a),sin(a),-sin(a),cos(a));	
}

vec2 iSphere( in vec3 ro, in vec3 rd, in vec4 sph )//from iq
{
	vec3 oc = ro - sph.xyz;
	float b = dot( oc, rd );
	float c = dot( oc, oc ) - sph.w*sph.w;
	float h = b*b - c;
	if( h<0.0 ) return vec2(-1.0);
	h = sqrt(h);
	return vec2(-b-h, -b+h );
}

float map(in vec3 p, vec2 sctime) {
	
	float res = 0.;
	
    vec3 c = p;
    c.xy = c.xy * sctime.x + vec2(c.y, c.x) * sctime.y;
	for (int i = 0; i < 10; ++i) {
        p =.7*abs(p)/dot(p,p) -.7;
        p.yz= csqr(p.yz);
        p=p.zxy;
        res += exp(-19. * abs(dot(p,c)));
        
	}
	return res/2.;
}



vec3 raymarch( in vec3 ro, vec3 rd, vec2 tminmax , vec2 sctime)
{
    tminmax += vec2(.5,.75) * sin( smoothTimeB * .3)*3.0;
   	vec3 one3 = vec3(1.,1.,1.);
    vec3 t = one3 * tminmax.x;
    
    vec3 dt = vec3(.07, 0.02, 0.05);
    vec3 col= vec3(0.);
    vec3 c = one3 * 0.;
    for( int i=0; i<60; i++ )
	{
     	vec3 s = vec3(2.0, 3.0, 4.0);   
        t+=dt*exp(-s*c);
        vec3 a = step(t,one3*tminmax.y);
        vec3 pos = ro+t*rd;
        
        c.x = map(ro+t.x*rd, sctime);
        c.y = map(ro+t.y*rd, sctime);
        c.z = map(ro+t.z*rd, sctime);               
        
        col = mix(col, .99*col+ .08*c*c*c, a);
    }
    
    vec3 c0 = vec3(0.4,0.3,0.99);
    vec3 c1 = vec3(0.9,0.7,0.0);
    vec3 c2 = vec3(0.9,0.1,0.2);
    return c0 * col.x + c1 * col.y + c2 * col.z;
}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

	float time = TIME;
    vec2 q = fragCoord.xy / RENDERSIZE.xy;
    vec2 p = -1.0 + 2.0 * q;
    p.x *= RENDERSIZE.x/RENDERSIZE.y;
    vec2 m = vec2(0.);
	if( iMouse.z>0.0 )m = iMouse.xy/RENDERSIZE.xy*3.14;
    m-=.5;

    // camera

    vec3 ro = zoom*vec3((inOut*sin(smoothTime*.125)+(1.385-zoomy))*3.);
    ro.yz*=rot(m.y+smoothTime*0.0725);
    ro.xz*=rot(m.x+ 0.01*smoothTime);
    vec3 ta = vec3( 0.0 , 0.0, 0.0 );
    vec3 ww = normalize( ta - ro );
    vec3 uu = normalize( cross(ww,vec3(0.0,1.0,0.0) ) );
    vec3 vv = normalize( cross(uu,ww));
    vec3 rd = normalize( p.x*uu + p.y*vv + 4.0*ww );

    
    vec2 tmm = iSphere( ro, rd, vec4(0.,0.,0.,2.) );

	// raymarch
    vec3 col = raymarch(ro,rd,tmm, vec2(sin(smoothTime*0.125), cos(smoothTimeC*0.25)));
    /*if (tmm.x<0.)col = texture(iChannel0, rd).rgb;
    else {
        vec3 nor=(ro+tmm.x*rd)/2.;
        nor = reflect(rd, nor);        
        float fre = pow(.5+ clamp(dot(nor,rd),0.0,1.0), 3. )*1.3;
        col += texture(iChannel0, nor).rgb * fre;
    
    }
	*/
	// shade
    
    col =  .5 *(log(1.+col));
    col = clamp(col,0.,1.);
    fragColor = vec4( col, 1.0 );

	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}