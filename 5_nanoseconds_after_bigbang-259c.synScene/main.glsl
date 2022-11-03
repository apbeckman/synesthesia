

// Created by Benoit Marini - 2020
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.


//Getting compil error since 11/2020 with 259 chars version, back to longer version
//-2 chars, thanks Xor !
float t2 = smoothTimeC;
vec4 renderMainImage() {
	vec4 o = vec4(0.0);
	vec2 F = _xy;

    vec2 R = RENDERSIZE.xy; 
    o-=o;
    for(float d,t = smoothTime*.1, i = 0. ; i > -1.; i -= .06 )          	
    {   d = fract( i -3.*t );                                          	
        vec4 c = vec4( ( F - R *.5 ) / R.y *d ,i,0 ) * 28.;            	
        for (int j=0 ; j++ <27; )                                      	
            c.xzyw = abs( c / dot(c,c)                                 	
                    -vec4( 7.-.2*sin(t) , 6.3 , .7 , 1.-cos(t/.8))/7.);	
       o*=pow(1.0+highhits) 
       o -= c * c.yzww  * d--*d  / vec4(3,5,1,1);                     
    }
    return o;
}

//Getting compil error since 11/2020 with 259 chars version, back to longer version
/*
vec4 renderMainImage() {
	vec4 o = vec4(0.0);
	vec2 F{ = _xy;

    vec2 R = RENDERSIZE.xy; 
    o-=o;
    for(float d,t = TIME*.1, i = 0. ; i > -1.; i -= .06 )          	
    {   d = fract( i -3.*t );                                          	
        vec4 c = vec4( ( F - R *.5 ) / R.y *d ,i,0 ) * 28.;            	
        for (int j=0 ; j++ <27; )                                      	
            c.xzyw = abs( c / dot(c,c)                                 	
                    -vec4( 7.-.2*sin(t) , 6.3 , .7 , 1.-cos(t/.8))/7.);	
       o += c * c.yzww  * (d-d*d)  / vec4(3,5,1,1);                     
    }
}*/
    
// 259 chars
// Thank you Xor !
/*
#define mainImage(o,F)                                                 	\
    vec2 R = RENDERSIZE.xy;                                           	\
    for(float d,t = TIME*.1, i = 0. ; i > -1.; i -= .06 )          	\
    {   d = fract( i -3.*t );                                          	\
        vec4 c = vec4( ( F - R *.5 ) / R.y *d ,i,0 ) * 28.;            	\
        for (int j=0 ; j++ <27; )                                      	\
            c.xzyw = abs( c / dot(c,c)                                 	\
                    -vec4( 7.-.2*sin(t) , 6.3 , .7 , 1.-cos(t/.8))/7.);	\
       o += c * c.yzww  * (d-d*d)  / vec4(3,5,1,1);                     \
    }
*/


// 262 chars
// Big thanks to FabriceNeyret2 for code reduction (and thus learning)
/*#define mainImage(o,F)                                               \
    vec2 R = RENDERSIZE.xy;                                           \
    for(float d, t = TIME*.1, i = 0. ; i > -1.; i -= .06 )            \
    {   d = fract( i -3.*t );                                          \
        vec4 c = vec4( ( F - R *.5 ) / R.y *d ,i,0 ) * 28.;            \
        for (int j=0 ; j++ <27; )                                      \
            c.xzyw = abs( c / dot(c,c)                                 \
                         - vec4( 1.-.03*sin(t) , .9 , .1 , .15 -.14*cos(t*1.3)) );\
       o += c * c.yzww * (d-d*d)  / vec4(3,5,1,1);                    \
}
/*

//original "short code" 298 chars
/*vec4 renderMainImage() {
	vec4 o = vec4(0.0);
	vec2 F{ = _xy;

    vec3 c;   
    float t = TIME*.1,i,d;    
	for(i=0.; i<1.; i+=.06)
    {
        d = fract(i+3.*t);   
        o = vec4( (F-RENDERSIZE.xy*.5)/RENDERSIZE.y*(1.-d) ,-i,0)*28.;   
    	for (int i=0 ; i++ <27;) o.xzyw = abs( o/dot(o,o) - vec4( 1.-.03*sin(t) , .9 , .1 , .15 -.14*cos(t*1.3)) );      
		c+= o.xyz*o.yzw*(d-d*d);
    }        
    o.rgb = c*vec3(.3,.2,1);
	return o; 
 } 

*/

// for more clear code see https://www.shadertoy.com/view/WtjyzR

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}