vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 


// Copyright Inigo Quilez, 2013 - https://iquilezles.org/
// I am the sole copyright owner of this Work.
// You cannot host, display, distribute or share this Work in any form,
// including physical and digital. You cannot use this Work in any
// commercial or non-commercial product, website or project. You cannot
// sell this Work and you cannot mint an NFTs of it.
// I share this Work for educational purposes, and you can link to it,
// through an URL, proper attribution and unmodified screenshot, as part
// of your educational material. If these conditions are too restrictive
// please contact me and we'll definitely work it out.
float smoothTime = (smooth_basstime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.125;
float smoothTimeB = (smooth_hightime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.125;
float smoothTimeC = (smooth_midtime * 0.75 + TIME * 0.125 + syn_Time * 0.25)*0.125;


float hash( float n )
{
    return fract(sin(n)*43758.5453);
}

float noise( in vec2 x )
{
    vec2 p = floor(x);
    vec2 f = fract(x);
    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y*57.0;
    return mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
               mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y);
}

vec2 map( vec2 p, in float offset )
{
	p.x += 0.1*sin( smoothTimeC + 2.0*p.y ) ;
	p.y += 0.1*sin( smoothTimeC + 2.0*p.x ) ;
	
	float a = noise(p*1.5 + sin(0.1*smoothTimeB))*6.2831;
	a -= offset;
	return vec2( cos(a), sin(a) );
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    vec2 p = fragCoord.xy / RENDERSIZE.xy;
	vec2 uv = -1.0 + 2.0*p;
	uv.x *= RENDERSIZE.x / RENDERSIZE.y;
		
    float offset = smoothTime + fragCoord.x/RENDERSIZE.x;
    
	float acc = 0.0;
	vec3  col = vec3(0.0);
	for( int i=0; i<64; i++ )
	{
		vec2 dir = map( uv, offset );
		
		float h = float(i)/64.0;
		float w = 4.0*h*(1.0-h);
		
		vec3 ttt = w*texture( image2, uv ).xyz;
		ttt *= mix( vec3(0.6,0.7,0.7), vec3(1.0,0.95,0.9), 0.5 - 0.5*dot( reflect(vec3(dir,0.0), vec3(1.0,0.0,0.0)).xy, vec2(0.707) ) );
		col += w*ttt;
		acc += w;
		
		uv += 0.006*dir;
	}
	col /= acc;
    
	float gg = dot( col, vec3(0.333) );
	vec3 nor = normalize( vec3( dFdx(gg), 0.5, dFdy(gg) ) );
	col += vec3(0.4)*dot( nor, vec3(0.7,0.01,0.7) );

	vec2 di = map( uv, offset );
	col *= 0.65 + 0.35*dot( di, vec2(0.707) );
	col *= 0.20 + 0.80*pow( 4.0*p.x*(1.0-p.x), 0.1 );
	col *= 1.7;

	fragColor = vec4( col, 1.0 );
	return fragColor; 
 } 



vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}