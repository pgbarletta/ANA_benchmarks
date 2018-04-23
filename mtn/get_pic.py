from pymol import cmd

cmd.color("skyblue", "1mtn")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_trace_mode",  1)
cmd.set("two_sided_lighting", "on")
cmd.set("reflect", 0)
cmd.set("ambient", 0.5)
cmd.set("ray_trace_mode",  0)
cmd.set('''ray_opaque_background''', '''off''')

cmd.png("pic_6_mtn_ch_differences.png", width=900, height=1100, dpi=600, ray=1)
