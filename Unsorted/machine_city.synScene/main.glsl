


// --- thanks to iq for all the ray-marching sdf info and awesomeness

vec3 op_u(vec3 d1, vec3 d2)
{
	return (d1.x < d2.x) ? d1 : d2;
}

void sphere_fold(inout vec3 p, inout float dr, float m_rad_sq, float f_rad_sq, float m_rad_sq_inv)
{
    float r_sq = dot(p, p);
    if (r_sq < m_rad_sq)
    {
        float t = f_rad_sq * m_rad_sq_inv;
        p *= t;
        dr *= t;
    }
    else if (r_sq < f_rad_sq)
    {
        float t = f_rad_sq / r_sq;
        p *= t;
        dr *= t;
    }
}

void box_fold(inout vec3 p, float fold_limit)
{
    p = clamp(p, -fold_limit, fold_limit) * 2.0 - p;
}

// estimators return (dist, mat_id, custom_value)

vec3 estimator_mandelbox(vec3 p, float scale, float m_rad_sq, float f_rad_sq, float fold_limit, float mat_id)
{
    vec3 off = p;
    float dr = 1.0+een;
    float dist = 1e20;
    float mrs_inv = 1.0 / m_rad_sq;
    for (int i = 0; i < 10; ++i)
    {
        box_fold(p, fold_limit);
        sphere_fold(p, dr, m_rad_sq, f_rad_sq, mrs_inv);

        p = scale * p + off;
        dr = dr * abs(scale) + 1.0;
        vec3 ot = p - vec3(0.5);
        dist = min(dist, dot(ot, ot));
    }
    return vec3(length(p) / abs(dr), mat_id, sqrt(dist));
}

vec3 mod_pos(vec3 p, float a, float b)
{
    p.zx = mod(p.zx, a) - b;  
    return p;
}

vec3 estimate(vec3 p)
{
    vec3 p_mb = mod_pos(p, 4.4, 2.2);
   	vec3 res_mb = estimator_mandelbox(p_mb, -2.5, 0.1, 2.5, 1.0, 0.2);
    // second
    vec3 p_pl = p;
    p_pl.y += 4.0;
    p_pl = mod_pos(p_pl, 2.0, 1.0);
    vec3 res_pl = estimator_mandelbox(p_pl, -1.5, 0.3, 2.9, 1.0, 0.1);

    return op_u(res_mb, res_pl);

}

vec4 ray_march(vec3 origin, vec3 direction)
{
    float total_distance = 0.0;
    int steps_total = 0;
    float orbit_trap;
    float material = -1.0;
    float dist;
    for (int steps = 0; steps < 100; ++steps)
    {
        vec3 point = origin + total_distance * direction;
        vec3 res = estimate(point);
        dist = res.x;
        
        if (dist < 0.003)
            break;
        
        total_distance += dist;
        steps_total++;
        
        orbit_trap = res.z;
        material = res.y;
    }
    float ao = float(steps_total) / 200.0;
    return vec4(total_distance, material, ao, orbit_trap);
}

vec3 compute_normal(vec3 pos)
{
    vec3 eps = vec3(0.001, 0.0, 0.0);
    vec3 pos00 = (pos + eps.xyy);
    vec3 pos01 = (pos - eps.xyy);
    vec3 pos10 = (pos + eps.yxy);
    vec3 pos11 = (pos - eps.yxy);
    vec3 pos20 = (pos + eps.yyx);
    vec3 pos21 = (pos - eps.yyx);
	vec3 normal = vec3(estimate(pos00).x - estimate(pos01).x, 
                       estimate(pos10).x - estimate(pos11).x, 
                       estimate(pos20).x - estimate(pos21).x);
	return normalize(normal);
}

vec3 render(vec3 origin, vec3 direction, vec2 uv)
{
    vec3 color = vec3(0.0, 0.0, 0.0);
    vec4 res = ray_march(origin, direction);
    
    vec3 light_dir = normalize(vec3(-0.6, -0.7, 0.5));
    float lg = length(uv - vec2(0.0));
    vec3 bg = exp(-vec3(5.0 - lg * 0.8, lg, lg * 0.5) * 2.0);
    color = bg;
    
    if (res.y > 0.0)
    {
        vec3 new_pos = origin + res.x * direction;
        float ao = 1.0 - res.z;   
        if (res.y == 0.2)
        {  
            vec3 normal = compute_normal(new_pos);
            
            vec3 h = normalize(light_dir - direction);
            float noh = clamp(dot(normal, h), 0.0, 1.0); 
            float voh = clamp(dot(-direction, h), 0.0, 1.0);
            float nol = clamp(dot(normal, light_dir), 0.0, 1.0);
            float nol_i = clamp(dot(normal, normalize(light_dir * vec3(-1.0, 0.0, -1.0))), 0.0, 1.0);

            float roughness = 0.54;
            float alphaSq = roughness * roughness * roughness * roughness;
            float denom = noh * noh * (alphaSq - 1.0) + 1.0;
            float D = alphaSq / (3.1416 * denom * denom);

            vec3 f0 = vec3(0.664, 0.824, 0.850);
            vec3 F = f0 + (1.0 - f0) * exp2((-5.55473 * voh - 6.98316) * voh);  

            vec3 spec = D * F;
            
            vec3 diff = nol * vec3(0.01, 0.49, 0.70) *  pow(vec3(ao), vec3(1.5, 1.2, 1.0));
            diff *= spec;
            diff += nol_i * vec3(0.41, 0.12, 0.16);
            diff *= vec3(cos(res.w * 0.71) * 0.5 + 0.5, sin(res.w * 0.01) * 0.01 + 0.5, 0.54);
            color = diff * ao;
        }
        else
        {
            vec3 diff = vec3(cos(res.w * 0.71) * 0.7 + 0.3, sin(res.w * 0.11) * 0.5 + 0.5, 0.54);
            color = diff * ao * ao * ao * vec3(1.64, 1.27, 0.99) * pow(vec3(ao),vec3(1.0, 1.2, 1.5));
        }
        
        float fog_factor = 0.1;
        float fog = 1.0 - exp(-res.x * fog_factor);
        color = mix(color, bg, fog);
    }
    
    return color;
}

vec4 renderMain() { 
 	vec4 out_FragColor = vec4(0.0);




    vec2 offsets[5];
    offsets[0] = vec2(0.0, 0.0);
    offsets[1] = vec2(0.25, 0.0);
    offsets[2] = vec2(-0.25, -0.0);
    offsets[3] = vec2(0.0, 0.25);
    offsets[4] = vec2(0.0, -0.25);
    
    vec2 mo = iMouse.xy / RENDERSIZE.xy;
    
    vec3 cam_pos = vec3(0.0  , sin(0.1  + 2.0 * mo.y), -smoothTime);
    vec3 cam_up = normalize(vec3(0.0, 0.8, 0.2));
    vec3 cam_right = vec3(1.0, 0.0, 0.0);
    vec3 cam_forward = normalize(cross(cam_up, cam_right));
    float focal_length = 2.77;
    vec3 color = vec3(0.0);
    
    const int numSamples = 1;
    for (int i = 0; i < numSamples; ++i)
    {
        vec2 coords_ss = (_xy.xy + offsets[i]) / RENDERSIZE.xy;
    	vec2 coords_cs = 2.0 * coords_ss - 1.0;
   	 	coords_cs.x *= RENDERSIZE.x / RENDERSIZE.y;
        
        vec3 ray_d = normalize(cam_forward * focal_length + cam_right * coords_cs.x + cam_up * coords_cs.y);
    	vec3 ray_o = cam_pos;
    
        color += render(ray_o, ray_d, coords_cs);
    }
    
    color /= float(numSamples);
	color = pow(color, vec3(0.4545));
    out_FragColor = vec4(color, 1.0);

return out_FragColor; 
 } 

