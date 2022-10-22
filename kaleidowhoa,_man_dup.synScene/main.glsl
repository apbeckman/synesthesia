float time = 0.0;
float rotTime = time_bass + syn_BeatTime*jump_rotate;
float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return min(max(d.x,max(d.y,d.z)),0.0) + length(max(d,0.0));
}

float opS( float d1, float d2 )
{
    return max(-d1,d2);
}

 mat3 rotationMatrix(vec3 axis, float angle)
{
    axis = normalize(axis);
    float s = sin(angle);
    float c = cos(angle);
    float oc = 1.0 - c;
    
    return mat3(oc * axis.x * axis.x + c,           oc * axis.x * axis.y - axis.z * s,  oc * axis.z * axis.x + axis.y * s,
                oc * axis.x * axis.y + axis.z * s,  oc * axis.y * axis.y + c,           oc * axis.y * axis.z - axis.x * s,
                oc * axis.z * axis.x - axis.y * s,  oc * axis.y * axis.z + axis.x * s,  oc * axis.z * axis.z + c);
}

float length8(vec3 x)
{
 	return pow(dot(vec3(1.0, 2.0, 1.0), pow(abs(x), vec3(8.0))), 1.0/8.0)  ;
}

float length8(vec2 x)
{
 	return pow(dot(vec2(1.0, 1.0), pow(abs(x), vec2(8.0))), 1.0/8.0)  ;
}

float sdTorus88( vec3 p, vec2 t )
{
    
  vec2 q = vec2(length8(p.xy)-t.x,p.z);
  return length8(q)-t.y;
}

float nsin(float x)
{
    return sin(x) * 0.5 + 0.5;
    
}

float map(vec3 p)
{

    vec3 q = p;

    float rep = 0.01;
        
    vec3 c = vec3(rep);
    p.z = mod(p.z,c.z)-0.5*c.z;

    
    vec3 p_s;
        
    float bars = 1000.0;
    float inner = 1000.0;
    float angle = 3.1415 * 0.33;
    
    float blockID = floor(q.z / rep);

    for ( int i = 0; i < 5; i ++)
    {
        
        p_s = p * rotationMatrix(vec3(0.0, 0.0, 1.0), angle * float(i));
        
        
        float cutout = 10000.;
        vec2 line = vec2( nsin(q.z*33.0), 0.005 * nsin(10.0 * q.z));
        
        p_s = p_s + vec3(
            sin(blockID * 11.0)* 0.1 + 0.3 -0.01*pow(syn_OnBeat,0.3)*beat_reaction,
            sin(q.z * sin(q.z+ rotTime* 0.01)) * sin(p.z* 4.0),
            0.0);
         
        p_s = p_s * rotationMatrix(vec3(0.0, 0.0, 1.0), 1.1 *  rotTime * 0.05  );
     	p_s = p_s * vec3(5.0*complexity, 1.0, 1.0);

        cutout = sdTorus88(p_s, line);
        
        inner = min(inner, cutout);

    }

        
    
    float result = inner;  
    return result;
}


void getCamPos(inout vec3 ro, inout vec3 rd)
{
    ro.z = time;
}

 vec3 gradient(vec3 p, float t) {
			vec2 e = vec2(0., t);

			return normalize( 
				vec3(
					map(p+e.yxx) - map(p-e.yxx),
					map(p+e.xyx) - map(p-e.xyx),
					map(p+e.xxy) - map(p-e.xxy)
				)
			);
		}


vec4 renderMainImage() {
	vec4 fragColor = vec4(0.0);
	vec2 fragCoord = _xy;

    if (_toPolarTrue(_uvc).y>1.0/6.0){
        discard;
    }
	time = time_forward*0.5;
    vec2 _p = (-RENDERSIZE.xy + 2.0*fragCoord.xy) / RENDERSIZE.y;
    vec3 ray = normalize(vec3(_p, 1.0));
    vec3 cam = vec3(0.0, 0.0, 0.0);
    bool hit = false;
    getCamPos(cam, ray);
    
    float depth = 0.1, d = 0.0, iter = 0.0, glow = 0.0;
    vec3 p;
    
    for( int i = 0; i < 80; i ++)
    {
    	p = depth * ray + cam;
        d = map(p);
                  
        if (d < 0.0005) {
			hit = true;
            break;
        }
        if ( depth > 20.0)
            break;      
        
        float ratio =  nsin(time * 0.1) * 0.01 + 0.03 + nsin(time)* 0.02;
		depth += d * ratio;
        glow += 0.01/sqrt(d);
		iter++;
                   
    }
    vec3 col = vec3(0.0);
    
    // if(hit)
    // 	col = vec3(1.0 - iter / 80.0);

    // // col = pow(col, vec3(
    // //     cos(floor(p.z *300.0 )) * 0.8 + 0.9 , 
    // //     0.9, 
    // //     sin(floor(p.z / 0.1)) * 0.3 + 0.4 ));
    
    fragColor = vec4(1.0 - iter / 80.0, depth, hit, glow);
    
	return fragColor; 
 } 

vec2 kaleidoscope(vec2 uvIn, float n) {
  vec2 uv = uvIn;
  float angle = PI/n;
  
  // uv = uv*vec2(RENDERSIZE.x/RENDERSIZE.y,1.0);
  float r = fract(length(uv));
  float a = atan(uv.y, uv.x)/angle;
  
  a = mix(fract(a), 1.0 - fract(a), mod(floor(a), 2.0))*angle;
  
  return vec2(cos(a), sin(a+0.05))*(r+0.03);
}

vec4 renderMain(){
	if(PASSINDEX == 0){
		return renderMainImage();
	}
    if(PASSINDEX == 1){
        // vec2 uv = _rotate(_uvc*0.5,TIME);
        vec2 uv = _rotate(_uvc, PI/2);
        float sides = 6.0 - 2.0*square + 2.0*octagon;
        uv = kaleidoscope(uv, sides);
        // uv = 
        uv *= 0.5;
        uv += 0.5;
        vec4 img = texture(backTex, uv);
        vec2 RENDERSIZE_FIRST = vec2(1408, 792);
        vec4 imgn = texture(backTex, uv+vec2(0,1)/RENDERSIZE_FIRST);
        vec4 imgs = texture(backTex, uv+vec2(0,-1)/RENDERSIZE_FIRST);
        vec4 imgw = texture(backTex, uv+vec2(1,0)/RENDERSIZE_FIRST);
        vec4 imge = texture(backTex, uv+vec2(-1,0)/RENDERSIZE_FIRST);


        vec3 finalCol = vec3(0.0);
        vec4 dat = (img*6.0+imgn*1.0+imgs*1.0+imgw*1.0+imge*1.0)/10.0;
        float glow = clamp(dat.a/12.0, 0.0, 1.0);
        float hit = dat.b;
        if (sunset>0.5){
            vec3 hsv = vec3(dat.g*1.0+time_color+color_selector+length(_uvc)*0.3+glow*0.25, 0.85, dat.r+glow);
            // finalCol.rgb = _rgb2hsv(finalCol.rgb);
            // finalCol.r = dat.r+time_color;
            // finalCol.r += length(_uvc)+finalCol.b;
            // finalCol.g += 0.5;
            finalCol = _hsv2rgb(hsv);
            finalCol.g -= pow((1.0-finalCol.r),1.5)*0.35+pow((1.0-finalCol.b), 1.5)*0.35;
            finalCol.g *= 2.0;
            if (finalCol.g > max(finalCol.b, finalCol.r)-glow*0.5){
                finalCol.g -= finalCol.g;
            }
            finalCol.b -= finalCol.r*0.4;
            finalCol.b *= 0.7;
            // finalCol.g += finalCol.b*0.2;
        }
        else if (two_color_minimal>0.5){
            finalCol = vec3(0.0);
            finalCol += (0.5+0.5*sin(dat.g*20.0+time_highs*4.0))*(1.0-hit)*col2;
            finalCol += glow*glow*1.5*col1;

        }
        else {
            finalCol.rgb = vec3(dat.r);
            finalCol.rgb = _rgb2hsv(finalCol.rgb);
            finalCol.r = glow*0.2+time_color*1.5+color_selector;
            finalCol.r += length(_uvc)+dat.g;
            finalCol.g += 0.8;
            finalCol.b += pow(glow,2.0);
            finalCol.b = clamp(finalCol.b, 0.0, 1.0);
            finalCol = _hsv2rgb(finalCol);
        }
        if (syn_MediaType > 0.5){
            vec3 imgCol = _loadUserImage(_rotate(vec2(1.0,1.0), 2*PI*(dat.r+dat.g + glow))*media_distort*media_distort).rgb;
            finalCol -= dot(imgCol, vec3(1.0))/3.0;
            finalCol = clamp(finalCol, 0.0, 1.0);
            finalCol += imgCol;
        }
        return vec4(finalCol, 1.0);
    }
}