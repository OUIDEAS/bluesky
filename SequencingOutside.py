import numpy as np
import matplotlib.pyplot as plt
import copy

# def step0(etas):
#     '''
#     CTAS Step 0: Verifying FCFS and position shifting if necessary
#     '''
#     gates = {'I': [], 'II': [], 'III': [], 'IV': []}

#     # Group aircraft by gate dynamically
#     for i in etas:
#         gates[i[3]].append(i)

#     # Print aircraft by gate
#     for gate, aircraft in gates.items():
#         print(f'GATE {gate} AIRCRAFT: {aircraft}')

#     # Sort lists dynamically and store in dictionaries
#     fcfs_sorted = {gate: sorted(aircraft, key=lambda x: x[2]) for gate, aircraft in gates.items()}
#     etas_sorted = {gate: sorted(aircraft, key=lambda x: x[0]) for gate, aircraft in gates.items()}


#     # for gate, aircraft in fcfs_sorted.items():
#     #     c=1
#     #     b=0
#         # for i in aircraft:
#         #     plt.scatter(b,c, marker = 's')
#         #     plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
#         #     c+=1
#         # plt.text(b, c+0.5, 'FCFS')
#         # b+=1
#         # c=1
#         # for i in etas_sorted[gate]:
#         #     plt.scatter(b,c, marker = 's')
#         #     plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
#         #     c+=1
#         # plt.text(b, c+0.5, 'Position Shift')
#         # plt.ylim((0, c+0.6))
#         # plt.xlim((-.1, 1.6))
#         # plt.title(f'Sorting for Gate{gate}')
#         # plt.show()
    
#     if fcfs_sorted == etas_sorted:
#         order = 'fcfs'
#         return fcfs_sorted, order
#     else:
#         order = 'pushback'
#         return etas_sorted, order
    
# def step12(fcfs, t_it, t_tran):
#     '''
#     Step 1: Applying in-trail constraints
#     '''
#     #In trail constraints
#     staff_it = {key: [vals[0]] if vals else [] for key, vals in fcfs.items()}

#     for key, vals in fcfs.items():
#         for i in range(1, len(vals)):
#             if staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
#                 staff_it[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0]])
#             else:
#                 staff_it[key].append([staff_it[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0]])
#     # print(fcfs)
#     # print(staff_it)
#     '''
#     Step 2: Applying transition constraints
#     '''
#     # rta = {'I': [], 'II': [], 'III': [], 'IV': []}
#     # ind = 0
#     # for key, vals in staff_it.items():
#     #     for i in vals:
#     #         rs = [i[0]+t_tran[ind][0], i[0]+t_tran[ind][1]]
#     #         print(t_tran[ind])
#     #         if ind <=1:
#     #             rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#     #         else:
#     #             rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+3}'])
#     #     ind+=1
#     rta = []
#     rta2 = []
#     ind = 0
#     # print(staff_it)
#     for key, vals in staff_it.items():
#         for i in vals:
#             print(vals)
#             print(i)
#             rs = [i[0]+i[4][0], i[0]+i[4][1]]
#             # print(t_tran[ind])
#         # if ind <=1:
#             rta.append([np.min(rs), i[1],i[2], i[3], f'R{np.argmin(rs)+1}', i[5]])
#             rta2.append([np.max(rs), i[1],i[2], i[3], f'R{np.argmax(rs)+1}', i[5]])
#             # else:
#             #     rta.append([np.min(rs), i[1],i[2], i[3], f'R{np.argmin(rs)+3}'])
#             #     rta2.append([np.max(rs), i[1],i[2], i[3], f'R{np.argmax(rs)+3}'])

#         ind+=1
#         # print(ind)

#             # if key == 'I' or key == 'II':
#             #     rs = [i[0]+t_tran[0], i[0]+t_tran[1]]
#             #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#             #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#             # if key == 'III' or key == 'IV':
#             #     rs = [i[0]+t_tran[2], i[0]+t_tran[3]]
#             #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#             #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#     cp = []
#     for i in rta:
#         cp.append(i)
#     cp = sorted(cp, key = lambda x: x[0])
#     cp = [[v[0], v[1], i, v[3], v[4], v[5]] for i, v in enumerate(cp)]
#     print(cp)
#     c=0
#     print(rta)


#     # for key, vals in staff_it.items():
        
#     #     for i in range(len(vals)):
#     #         if c<=4:
#     #             plt.scatter(c, fcfs[key][i][0], marker = 's')
#     #             plt.scatter(c+1, vals[i][0], marker = 's')
                
#     #             plt.text(c-0.1, fcfs[key][i][0], f'{fcfs[key][i][1]}')
#     #             plt.text(c+1-0.1, vals[i][0], f'{vals[i][1]}')

#     #             plt.plot([c, c+1], [fcfs[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
#     #         # plt.plot([c+1, c+2], [vals[i][0], rta[key][i][0]], linestyle = '--', color = 'k')
#     #         else:
#     #             plt.scatter(c, fcfs[key][i][0], marker = 's')
#     #             plt.scatter(c-1, vals[i][0], marker = 's')
                
#     #             plt.text(c+0.1, fcfs[key][i][0], f'{fcfs[key][i][1]}')
#     #             plt.text(c-1-0.1, vals[i][0], f'{vals[i][1]}')

#     #             plt.plot([c, c-1], [fcfs[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
#     #     c+=2
#     #     if c == 4:
#     #         c+=4

#     # c=0
#     # for i in cp:
#     #     print(f'{i[2]} ', end='')
#     #     plt.scatter(5, i[0], marker = 's')
#     #     # if 'A' in i[2]:
#     #         # plt.text(0.9, i[0], f'{i[2]}', weight = 'bold')
#     #     # else:
#     #     plt.text(5, i[0], f'{i[2]}', weight = 'bold')
#     #     c+=1
#     # # flat = np.array(rta).flatten()
#     # plt.show()


#     rta_p = {'R1': [], 'R2': []}
#     for i in rta:
#         rta_p[i[4]].append(i)
#     print(rta_p)

#     non_rta = {'R1': [], 'R2': []}
#     for i in rta2:
#         non_rta[i[4]].append(i)
#     print(non_rta)
#             # print(ind)

#                 # if key == 'I' or key == 'II':
#                 #     rs = [i[0]+t_tran[0], i[0]+t_tran[1]]
#                 #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#                 #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#                 # if key == 'III' or key == 'IV':
#                 #     rs = [i[0]+t_tran[2], i[0]+t_tran[3]]
#                 #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
#                 #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
#         # print(rta)
#     return cp, rta_p, non_rta

# def step3(cp, rta_p, non_rta, order, t_it):
#     '''
#     Step 3: Calculating STA
#     '''

#     sta_p = {key: [vals[0]] if vals else [] for key, vals in rta_p.items()}
#     non_sta = {key: [vals[0]] if vals else [] for key, vals in non_rta.items()}

#     #staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
#     # print(sta_p.keys())
#     for key, vals in rta_p.items():
#         for i in range(1, len(vals)):
#             if sta_p[key][i-1] and vals[i][0] > sta_p[key][i-1][0]+t_it:
#                 sta_p[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5]])
#             else:
#                 sta_p[key].append([sta_p[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5]])
#     # print(sta_p)

#     for key, vals in non_rta.items():
#         for i in range(1, len(vals)):
#             if non_sta[key][i-1] and vals[i][0] > non_sta[key][i-1][0]+t_it:
#                 non_sta[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5]])
#             else:
#                 non_sta[key].append([non_sta[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5]])

#     c=0
#     for key, vals in rta_p.items():
#         for i in vals:
#                 print(i)
#                 plt.scatter(c, i[0], marker = 's')
#                 plt.text(c-0.1, i[0], f'{i[1]}')
#         plt.text(c-0.2, i[0]+15, f'RTA Runway {key}')
#         c+=2
#     plt.show()

#     return sta_p, non_sta

# def step4(sta_p, dtmax):
#     '''
#     CTAS Step 4 Delay Distribution Function
#     '''
#     delay = []
#     ddfC = []
#     ddfT = []

#     for key, vals in sta_p.items():
#         for i in vals:
#             delay.append(i[0]-i[5], i[0], i[1], i[2], i[3], i[4], i[5])

#     for i in delay:
#         if i[0] <= dtmax:
#             ddfC.append([0, i[1], i[2], i[3], i[4], i[5], i[6]])
#             ddfT.append(i)
#         else:
#             ddfC.append([i[0]-dtmax, i[1], i[2], i[3], i[4], i[5], i[6]])
#             ddfT.append([dtmax, i[1], i[2], i[3], i[4], i[5], i[6]])


#     return delay, ddfC, ddfT

# def step5():

#     return

def sort(etas, t_it, t_tran, dtmax):

    '''
    CTAS Step 0: Verifying FCFS and position shifting if necessary
    '''
    gates = {'I': [], 'II': [], 'III': [], 'IV': []}

    # Group aircraft by gate dynamically
    for i in etas:
        gates[i[3]].append(i)

    # Print aircraft by gate
    for gate, aircraft in gates.items():
        print(f'GATE {gate} AIRCRAFT: {aircraft}')

    # Sort lists dynamically and store in dictionaries
    fcfs_sorted = {gate: sorted(aircraft, key=lambda x: x[2]) for gate, aircraft in gates.items()}
    etas_sorted = {gate: sorted(aircraft, key=lambda x: x[0]) for gate, aircraft in gates.items()}


    for gate, aircraft in fcfs_sorted.items():
        c=1
        b=0
        for i in aircraft:
            plt.scatter(b,c, marker = 's')
            plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
            c+=1
        plt.text(b, c+0.5, 'FCFS')
        b+=1
        c=1
        for i in etas_sorted[gate]:
            plt.scatter(b,c, marker = 's')
            plt.text(b+0.005, c, f'{i[1]}, {i[0]:.2f}, {i[2]}', weight = 'bold')
            c+=1
        plt.text(b, c+0.5, 'Position Shift')
        plt.ylim((0, c+0.6))
        plt.xlim((-.1, 1.6))
        plt.title(f'Sorting for Gate {gate}')
        plt.show()
    
    if fcfs_sorted == etas_sorted:
        order = 'fcfs'
        # return fcfs_sorted, order
    else:
        order = 'pushback'
        # return etas_sorted, order

    '''
    Step 1: Applying in-trail constraints
    '''
    #In trail constraints
    if order == 'fcfs':
        print('FCFS SORTED:', fcfs_sorted,'\n\n')
        staff_it = {key: [vals[0]] if vals else [] for key, vals in fcfs_sorted.items()}

        for key, vals in fcfs_sorted.items():
            for i in range(1, len(vals)):
                if staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
                    staff_it[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0], vals[i][-1]])
                else:
                    staff_it[key].append([staff_it[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0], vals[i][-1]])
    
    else:
        print('ETA SORTED:', etas_sorted,'\n\n')
        staff_it = {key: [vals[0]] if vals else [] for key, vals in etas_sorted.items()}

        for key, vals in etas_sorted.items():
            for i in range(1, len(vals)):
                if staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
                    staff_it[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0], vals[i][-1]])
                else:
                    staff_it[key].append([staff_it[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][0], vals[i][-1]])

    '''
    Step 2: Applying transition constraints
    '''
    print('STAFF_IT:', staff_it, '\n\n')
    rta = []
    rta2 = []
    for key, vals in staff_it.items():
        for i in vals:
            rs = [i[-1][0]+i[4][0], i[-1][1]+i[4][1]]
            # retas = [i[-1][0], i[-1][1]]
            # print(retas)
            # print(rs)
            # print(retas[np.argmin(rs)], retas[rs.index(np.min(rs))])
            # print(retas[rs.index(np.min(rs))]+np.min(rs), f'R{np.argmin(rs)+1}')
            if np.argmin(rs)+1 == 2 and i[3] == 'I':
                i[3] = 'II'
            elif np.argmin(rs)+1 == 2 and i[3] == 'IV':
                i[3] = 'III'
            elif np.argmin(rs)+1 == 1 and i[3] == 'II':
                i[3] = 'I'
            elif np.argmin(rs)+1 == 1 and i[3] == 'III':
                i[3] = 'IV'
            rta.append([np.min(rs), i[1],i[2], i[3], f'R{np.argmin(rs)+1}', i[5], i[4], i[0], i[6]])
            rta2.append([np.max(rs), i[1],i[2], i[3], f'R{np.argmax(rs)+1}', i[5], i[4], i[0], i[6]])

    cp = []
    for i in rta:
        cp.append(i)
    cp2 = sorted(cp, key = lambda x: x[0])
    cp3 = [[v[0], v[1], i, v[3], v[4], v[5], v[6], v[7], v[8]] for i, v in enumerate(cp2)]
    print('CP:',cp3, '\n\n')
    print('RTA:', rta, '\n\n')

    # for key, vals in staff_it.items():
        
    #     for i in range(len(vals)):
    #         if c<=4:
    #             plt.scatter(c, fcfs[key][i][0], marker = 's')
    #             plt.scatter(c+1, vals[i][0], marker = 's')
                
    #             plt.text(c-0.1, fcfs[key][i][0], f'{fcfs[key][i][1]}')
    #             plt.text(c+1-0.1, vals[i][0], f'{vals[i][1]}')

    #             plt.plot([c, c+1], [fcfs[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
    #         # plt.plot([c+1, c+2], [vals[i][0], rta[key][i][0]], linestyle = '--', color = 'k')
    #         else:
    #             plt.scatter(c, fcfs[key][i][0], marker = 's')
    #             plt.scatter(c-1, vals[i][0], marker = 's')
                
    #             plt.text(c+0.1, fcfs[key][i][0], f'{fcfs[key][i][1]}')
    #             plt.text(c-1-0.1, vals[i][0], f'{vals[i][1]}')

    #             plt.plot([c, c-1], [fcfs[key][i][0], vals[i][0]], linestyle = '--', color = 'k')
    #     c+=2
    #     if c == 4:
    #         c+=4

    # c=0
    # for i in cp:
    #     print(f'{i[2]} ', end='')
    #     plt.scatter(5, i[0], marker = 's')
    #     # if 'A' in i[2]:
    #         # plt.text(0.9, i[0], f'{i[2]}', weight = 'bold')
    #     # else:
    #     plt.text(5, i[0], f'{i[2]}', weight = 'bold')
    #     c+=1
    # # flat = np.array(rta).flatten()
    # plt.show()


    rta_p = {'R1': [], 'R2': []}
    for i in cp3:
        rta_p[i[4]].append(i)

    print('RTA_P:\n\n',rta_p, '\n\n')

    non_rta = {'R1': [], 'R2': []}
    for i in rta2:
        non_rta[i[4]].append(i)
    # print(non_rta)
            # print(ind)

                # if key == 'I' or key == 'II':
                #     rs = [i[0]+t_tran[0], i[0]+t_tran[1]]
                #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
                #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
                # if key == 'III' or key == 'IV':
                #     rs = [i[0]+t_tran[2], i[0]+t_tran[3]]
                #     # print(np.min(rs), np.argmin(rs), rs.index(np.min(rs)))
                #     rta[key].append([np.min(rs), i[2],i[1], i[3], f'R{np.argmin(rs)+1}'])
        # print(rta)

    '''
    Step 3: Calculating STA
    '''

    sta_p = {key: [vals[0]] if vals else [] for key, vals in rta_p.items()}
    print('INITIAL STA', sta_p, '\n\n')
    non_sta = {key: [vals[0]] if vals else [] for key, vals in non_rta.items()}

    #staff_it[key][i-1] and vals[i][0]>= staff_it[key][i-1][0]+t_it:
    # print(sta_p.keys())
    for key, vals in rta_p.items():
        for i in range(1, len(vals)):
            if sta_p[key][i-1] and vals[i][0] > sta_p[key][i-1][0]+t_it:
                sta_p[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], vals[i][8]])
            else:
                sta_p[key].append([sta_p[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], vals[i][8]])
    print('STA_P', sta_p, '\n\n')

    for key, vals in non_rta.items():
        for i in range(1, len(vals)):
            if non_sta[key][i-1] and vals[i][0] > non_sta[key][i-1][0]+t_it:
                non_sta[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], vals[i][8]])
            else:
                non_sta[key].append([non_sta[key][i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], vals[i][8]])

    c=0
    for key, vals in rta_p.items():
        for i in vals:
                print(i)
                plt.scatter(c, i[0], marker = 's')
                plt.text(c-0.1, i[0], f'{i[1]}')
        if key[0]:
            plt.text(c-0.2, i[0]+15, f'RTA Runway {key}')
        c+=2
    plt.show()

    '''
    CTAS Step 4 Delay Distribution Function
    '''
    delay = {'R1':[], 'R2':[]}
    ddfC = {'R1':[], 'R2':[]}
    ddfT = {'R1':[], 'R2':[]}

    for key, vals in sta_p.items():
        for i in vals:
            delay[key].append([i[0]-i[5], i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8]])

    for key, vals in delay.items():
        for i in vals:
            if i[0] <= dtmax:
                ddfC[key].append([0, i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]])
                ddfT[key].append(i)
            else:
                ddfC[key].append([i[0]-dtmax, i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]])
                ddfT[key].append([dtmax, i[1], i[2], i[3], i[4], i[5], i[6], i[7], i[8], i[9]])

    '''
    CTAS Step 5 Scheduled Time of Arrival at Meter Gates
    '''

    staff = {'R1': [], 'R2': []}
    
    for key, vals in sta_p.items():
        for i in range(len(vals)):
            # print(c)
            # print(vals, vals[i])
            rs = [vals[i][6][0], vals[i][6][1]]
            rs_p = np.min(rs)
            # print(rs_p, rs)
            if i == 0:
                staff[key].append([vals[i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], ddfC[key][i][0], vals[i][-1]])
            # elif i == 1:
            #     staff[key].append([vals[i][7]+ddfC[key][i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], ddfC[key][i][0], vals[i][-1]])
            else:
                if vals[i-1][0]+t_it > vals[i][7]+ddfC[key][i][0]:
                    staff[key].append([vals[i-1][0]+t_it, vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], ddfC[key][i][0], vals[i][-1]])
                else:
                    staff[key].append([vals[i][7]+ddfC[key][i][0], vals[i][1], vals[i][2], vals[i][3], vals[i][4], vals[i][5], vals[i][6], vals[i][7], ddfC[key][i][0], vals[i][-1]])
                if staff[key][i][0] - staff[key][i-1][0] < t_it:
                    staff[key][i][0] = staff[key][i][0] + (t_it - (staff[key][i][0] - staff[key][i-1][0]))

    print('STAFF_IT:', staff_it, '\n\n')
    print('STA:', sta_p, '\n\n')
    print('STAFF:', staff, '\n\n')
    staff2 = copy.deepcopy(staff)
    
    run = []
    j = 0
    c=0
    while c == 0:
        try:
            if staff['R1'][0][0] < staff['R2'][0][0]:
                # print(staff['R1'][0], '\n')
                run.append(staff['R1'][0])
                # print(run, '\n')
                staff['R1'].pop(0)
                # print(staff['R1'], '\n')
                j = 'R1'
            else:
                # print(staff['R2'][0], '\n')
                run.append(staff['R2'][0])
                # print(run, '\n')
                staff['R2'].pop(0)
                # print(staff['R2'], '\n')
                j = 'R2'
        except:
            # if j == 'R1':
            #     run.append(staff['R2'][0])
            #     staff['R2'].pop(0)
            # else:
            #     run.append(staff['R1'][0])
            # c=1
            if staff['R1']:
                for i in range(len(staff['R1'])):
                    run.append(staff['R1'][0])
                    staff['R1'].pop(0)
            elif staff['R2']:
                for i in range(len(staff['R2'])):
                    run.append(staff['R2'][0])
                    staff['R2'].pop(0)
            c=1
    c = 0
    # if run[0][0] <run[1][0] and run[1][0] - run[0][0] < t_it:
    #     run[1][0]+=t_it
        
    for i in cp3:
        # print(f'{i[2]} ', end='')
        plt.scatter(1, i[0], marker = 's')
        plt.text(0.995, i[0], f'{i[1]}', weight = 'bold')
        plt.text(1, i[0], f'{i[4]}', weight = 'bold')
        c+=1
    plt.text(1, cp3[0][0] - 100, 'RTA')
    # print('\n\nRUN LIST',run)
    print(staff)
    for i in run:
        # print(f'{i[2]} ', end='')
        plt.scatter(1.05, i[0], marker = 's')
        plt.text(1.045, i[0], f'{i[1]}', weight = 'bold')
        plt.text(1.05, i[0], f'{i[4]}', weight = 'bold')
        c+=1
    plt.text(1.05, cp3[0][0]-100, 'Step 5\nRTA' )
    for key, vals in sta_p.items():
        for i in vals:
            plt.scatter(1.025, i[0], marker = 's')
            plt.text(1.02, i[0], f'{i[1]}', weight = 'bold')
            plt.text(1.025, i[0], f'{i[4]}', weight = 'bold')

    plt.text(1.02, cp3[0][0]-100, 'Calculated\nSTA')
                # plt.text(c-0.1, i[0], f'{i[1]}')
        # plt.text(c-0.2, i[0]+15, f'STA Runway {key}')
    plt.xlim([0.95, 1.15])
    plt.show()

    return sta_p, run, cp3, staff2