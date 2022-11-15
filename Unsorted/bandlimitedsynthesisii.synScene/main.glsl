//vec4 iMouse = vec4(MouseXY*RENDERSIZE, MouseClick, MouseClick); 

			//******** Common Code Begins ********

// https://www.shadertoy.com/view/wsjfRD
// A white noise function.
float rand(vec3 p) {
    return fract(sin(dot(p, vec3(12.345, 67.89, 412.12))) * 42123.45) * 2.0 - 1.0;
}

float perlin(vec3 p) {
    vec3 u = floor(p);
    vec3 v = fract(p);
    vec3 s = smoothstep(0.0, 1.0, v);
    
    float a = rand(u);
    float b = rand(u + vec3(1.0, 0.0, 0.0));
    float c = rand(u + vec3(0.0, 1.0, 0.0));
    float d = rand(u + vec3(1.0, 1.0, 0.0));
    float e = rand(u + vec3(0.0, 0.0, 1.0));
    float f = rand(u + vec3(1.0, 0.0, 1.0));
    float g = rand(u + vec3(0.0, 1.0, 1.0));
    float h = rand(u + vec3(1.0, 1.0, 1.0));
    
    return mix(mix(mix(a, b, s.x), mix(c, d, s.x), s.y),
               mix(mix(e, f, s.x), mix(g, h, s.x), s.y),
               s.z);
}


// The MIT License
// Copyright © 2020 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// A simple way to prevent aliasing of cosine functions (the color
// palette in this case is made of 8 layers) by attenuating them
// when their oscillations become smaller than a pixel. Left is
// direct use of cos(x), right is band-limited cos(x).
//
// Box-filtering of cos(x):
//
// (1/w)∫cos(t)dt with t ∈ (x-½w, x+½w)
// = [sin(x+½w) - sin(x-½w)]/w
// = cos(x)·sin(½w)/(½w)
//
// Can approximate smoothstep(2π,0,w) ≈ sin(w/2)/(w/2),
// which you can also see as attenuating cos(x) when it 
// oscilates more than once per pixel. More info:
//
// https://iquilezles.org/www/articles/bandlimiting/bandlimiting.htm
//
// Related Shader:
//   https://www.shadertoy.com/view/WtScDt
//   https://www.shadertoy.com/view/wtXfRH
//   https://www.shadertoy.com/view/3tScWd

bool mode;

vec3 fcos( vec3 x )
{
    if( mode) return cos(x);                // naive

    vec3 w = fwidth(x);
    #if 0
    return cos(x) * sin(0.5*w)/(0.5*w);     // filtered-exact
	#else
    return cos(x) * smoothstep(6.28,0.0,w); // filtered-approx
	#endif  
}

vec3 getColor( in float t )
{
    vec3 col = vec3(0.4,0.4,0.4);
    col += 0.12*fcos(6.28318*t*  1.0+vec3(0.0,0.8,1.1));
    col += 0.11*fcos(6.28318*t*  3.1+vec3(0.3,0.4,0.1));
    col += 0.10*fcos(6.28318*t*  5.1+vec3(0.1,0.7,1.1));
    col += 0.09*fcos(6.28318*t*  9.1+vec3(0.2,0.8,1.4));
    col += 0.08*fcos(6.28318*t* 17.1+vec3(0.2,0.6,0.7));
    col += 0.07*fcos(6.28318*t* 31.1+vec3(0.1,0.6,0.7));
    col += 0.06*fcos(6.28318*t* 65.1+vec3(0.0,0.5,0.8));
    col += 0.06*fcos(6.28318*t*115.1+vec3(0.1,0.4,0.7));
    col += 0.09*fcos(6.28318*t*265.1+vec3(1.1,1.4,2.7));
    return col;
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    // coordiantes
	vec2 p = (2.0*fragCoord-RENDERSIZE.xy)/RENDERSIZE.y;
    vec2 w = p;
    
    // separation
    //float th = (iMouse.z>0.001) ? (2.0*iMouse.x-RENDERSIZE.x)/RENDERSIZE.y : 1.8*sin(TIME);
    //mode = (w.x-th<0.0);
    
    // deform 1
    p *=( 0.125); // bigger number = smaller center circle
    p.xy+=_uvc*FOV*PI;
    p = 0.5*p/dot(p,p);
        p.xy *=1.0+((p.xy*Melt*PI*_uvc)/_xy);
    p +=_rotate(p*_uvc*PI, Rotate)*Rotate;

    vec2 q = p;
    p.x += ((smoothTime) * (rate*0.5));
    
    // deform 2
    p += (0.2+(warp1 + sin(smoothTimeC * 0.061) ) )*cos( 1.5*p.yx +( 0.06*1.0*(sin(smoothTimeB*0.2) +  bassBoost))  + vec2(0.1,1.1) );
	p += (0.2+(warp2 + sin(smoothTimeB * 0.1)))*cos( 2.4*p.yx + 0.06*1.6*(smoothTimeC) + vec2(4.5,2.6) );
	p += (0.2+warp3)*cos( 3.3*p.yx + 0.06*1.2*smoothTimeB + vec2(3.2,3.4) );
	p += (0.2+warp3)*cos( 4.2*p.yx + 0.06*1.7*(smoothTimeC) + vec2(1.8,5.2) );
	p += 0.2*cos( 9.1*p.yx + 0.06*1.1*smoothTimeB + vec2(6.3,3.9) );
    // base color pattern
    vec3 col = getColor( 0.5*length(p) );
    
    // lighting
    col *= 1.4 - 0.07*length(q);

    // separation
   // col *= smoothstep(0.005,0.010,abs(w.x-th));
    
    // palette
    if( w.y<-1.0 ) col = getColor( fragCoord.x/RENDERSIZE.x );
 
   fragColor = vec4( col, 1.0 );
	return fragColor; 
 } 


vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}