// String Theory by nimitz (twitter: @stormoid)
// https://www.shadertoy.com/view/XdSSz1
// License Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
// Contact the author for other licensing options

#define BASE_ANGLE PI
#define ANGLE_DELTA 0.012
#define XOFF 1.7

#define time TIME*0.5
mat2 mm2(in float a){float c = cos(a), s = sin(a);return mat2(c,-s,s,c);}


float f(vec2 p, float featureSize)
{
	p.x = sin(p.x*1.+((bass_time)+((bass_time))))*sin(p.x*0.1+(bass_time*0.5))*2.1;	
    p += sin(p.x*(1.5+complexity)*(.1+complexity*0.5));
    return smoothstep(-0.0,featureSize*(1.0+syn_Intensity*syn_MidHighLevel),abs(p.y));
}

vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    float aspect = RENDERSIZE.x/RENDERSIZE.y;
    float featureSize = 60./((RENDERSIZE.x*aspect+RENDERSIZE.y));

    vec2 p = fragCoord.xy / RENDERSIZE.xy*(6.5+zoom)-(3.25+zoom*0.5);
	p.x *= aspect;
	p.y = abs(p.y+_uvc.y);

	vec3 col = vec3(0);
	for(float i=0.;i<16.;i++)
	{
		vec3 col2 = (sin(vec3(3.75,2.65,2.2)*(1.0+syn_Intensity*syn_HighLevel)+i*0.25)*0.5+0.54)*(1.-f(p,featureSize))*(1.0+syn_Intensity*syn_HighLevel);
		col = max(col,col2);
		
        p.x -= (XOFF+expand);
        p.y -= 1.5+yPos;
		p*= mm2(i*(angle+ANGLE_DELTA)+(BASE_ANGLE+angle2));
		
        vec2 pa = vec2(abs(p.x-.9),abs(p.y));
        vec2 pb = vec2(p.x,abs(p.y));
        
        p = mix(pa,pb,smoothstep(-.07,(.07+(0.01*syn_BassLevel)),sin(smoothTime*0.0024)+.1));
	}
	fragColor = vec4(col,1.0);
    return fragColor;
}
vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
}