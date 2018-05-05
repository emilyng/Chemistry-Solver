
import sympy

HI, HII, HM, HeI, HeII, HeIII, H2I, H2II, de = sympy.sympify("HI, HII, HM, HeI, HeII, HeIII, H2I, H2II, de")
k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, k32 = sympy.sympify("k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, k25, k26, k27, k28, k29, k30, k31, k32")    
  
r1 = (HI + de), (HII + de + de), k1
r2 = (HII + de), (HI), k2
r3 = (HeI + de), (HeII + de + de), k3
r4 = (HeII + de), (HeI), k4
r5 = (HeII + de), (HeIII + de + de), k5
r6 = (HeIII + de), (HeII), k6
r7 = (HI + de), (HM), k7
r8 = (HM + HI), (H2I + de), k8
r9 = (HI + HII), (H2II), k9
r10 = (H2II + HI), (H2I + HII), k10
r11 = (H2I + HII), (H2II + HI), k11
r12 = (H2I + de), (2*HI + de), k12
r13 = (H2I + HI), (3*HI), k13
r14 = (HM + de), (HI + 2*de), k14
r15 = (HM + HI), (2*HI + de), k15
r16 = (HM + HII), (2*HI), k16
r17 = (HM + HII), (H2II + de), k17
r18 = (H2II + de), (2*HI), k18
r19 = (H2II + HM), (HI + H2I), k19

all_reactions = [r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11 ,r12, r13, r14, r15, r16, r17, r18, r19]


'''
r7 = (HI + HI), (HII + de + HI), k7
r8 = (HI + HeI), (HII + de + HeI), k8    
r9 = (HI), (HII + de), k9
r10 = (HeI), (HeII + de), k10
r11 = (HeII), (HeIII + de), k11
r12 = (HI + de), (HM), k12
r13 = (HM + HI), (H2I + de), k13
r14 = (HI + HII), (H2II + de), k14
r15 = (H2II + HI), (H2I + HII), k15
r16 = (H2I + HII), (H2II + HI), k16
r17 = (H2I + de), (HI + HI + de), k17
r18 = (H2I + HI), (HI + HI + HI), k18
r19 = (HM + de), (HI + de + de), k19
r20 = (HM + HI), (HI + HI + de), k20
r21 = (HM + HII), (HI + HI), k21
r22 = (HM + HII), (H2II + de), k22
r23 = (H2II + de), (HI + HI), k23
r24 = (H2II + HM), (H2I + HI), k24
r25 = (HI + HI + HI), (H2I + HI), k25
r26 = (HI + HI + H2I), (H2I + H2I), k26
r27 = (HM), (HI + de), k27
r28 = (H2II), (HI + HII), k28
r29 = (H2I), (H2II + de), k29
r30 = (H2II), (HI + HI), k30
r31 = (H2I), (HI + HI), k31
r32 = (HI + HI), (H2I), k32 #grain?
'''