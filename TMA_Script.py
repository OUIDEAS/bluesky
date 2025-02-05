# %%
import numpy as np
import matplotlib.pyplot as plt
import copy
etas = []
for i in range(0, 5):
    time = np.random.uniform(low = 12, high = 30)
    fcfs = np.random.uniform(low = 1, high = 6)
    etas.append([int(fcfs), time, i, f'A{i}'])
print(etas)
fcfs_etas = sorted([i[:] for i in etas], key=lambda x: x[0])

c = 1
for i in fcfs_etas:
    i.append(c)
    plt.scatter(0, c)
    plt.text(0.005, c, f'AC{i[2]}, {i[1]:.2f}, {i[0]}', weight = 'bold')
    c+=1
plt.text(0, 5.5, 'FCFS', weight = 'bold')
print('FCFS:', fcfs_etas)

eta_sort_a = sorted([i[:] for i in etas], key=lambda x: x[1])
c=1
# new_eta = [i[:] for i in eta_sort]
for i in eta_sort_a:
    i.append(c)
    plt.scatter(1, c)
    plt.text(1.005, c, f'AC{i[2]}, {i[1]:.2f}, {i[0]}', weight = 'bold')
    c+=1
plt.text(1, 5.5, 'Position Shift', weight = 'bold')

index_eta = sorted(eta_sort_a, key = lambda x: x[2])

index_fcfs = sorted(fcfs_etas, key = lambda x: x[2])

print('ETA Sorted:', eta_sort_a)

for i in range(0, len(index_fcfs)):
    plt.plot([0, 1], [index_fcfs[i][4], index_eta[i][4]], color = 'black', linestyle = '--')
plt.xlim([-.25, 1.5])
plt.show()



# %%
etas_b = []
for i in range(0, 5):
    time = np.random.uniform(low = 12, high = 30)
    fcfs = np.random.uniform(low = 1, high = 6)
    etas_b.append([int(fcfs), time, i, f'B{i}'])

fcfs_etas_b = sorted([i[:] for i in etas_b], key=lambda x: x[0])

c = 1
for i in fcfs_etas_b:
    i.append(c)
    plt.scatter(0, c)
    plt.text(0.005, c, f'AC{i[2]}, {i[1]:.2f}, {i[0]}', weight = 'bold')
    c+=1
plt.text(0, 5.5, 'FCFS', weight = 'bold')
print('FCFS:', fcfs_etas_b)

eta_sort_b = sorted([i[:] for i in etas_b], key=lambda x: x[1])
c=1
# new_eta = [i[:] for i in eta_sort]
for i in eta_sort_b:
    i.append(c)
    plt.scatter(1, c)
    plt.text(1.005, c, f'AC{i[2]}, {i[1]:.2f}, {i[0]}', weight = 'bold')
    c+=1
plt.text(1, 5.5, 'Position Shift', weight = 'bold')
print('ETA Sorted:', eta_sort_b)

index_eta_b = sorted(eta_sort_b, key = lambda x: x[2])

index_fcfs_b = sorted(fcfs_etas_b, key = lambda x: x[2])

for i in range(0, len(index_fcfs_b)):
    plt.plot([0, 1], [index_fcfs_b[i][4], index_eta_b[i][4]], color = 'black', linestyle = '--')
plt.xlim([-.25, 1.5])
plt.show()

# %%
'''
CTAS Steps 1 & 2 In-Trail Constrains and Runway Threshold Landing Order
'''
in_trail = 5
transition_time = 30
#i[1] = eta
staff_it_a = []
staff_it_a.append([eta_sort_a[0][1], eta_sort_a[0][2], eta_sort_a[0][3]])
for i in range(1, len(eta_sort_a)):
    if eta_sort_a[i][1] >= staff_it_a[i-1][0]+in_trail:
        staff_it_a.append([eta_sort_a[i][1], eta_sort_a[i][2], eta_sort_a[i][3]])
    else:
        staff_it_a.append([staff_it_a[i-1][0]+in_trail, eta_sort_a[i][2], eta_sort_a[i][3]])
# print(eta_sort_a)
print('STAFF_it A:', staff_it_a)
rta_a = [[i[0]+transition_time, i[1], i[2], i[0]] for i in staff_it_a]
print('RTA A:', rta_a)

staff_it_b = []
staff_it_b.append([eta_sort_b[0][1], eta_sort_b[0][2], eta_sort_b[0][3]])
for i in range(1, len(eta_sort_b)):
    if eta_sort_b[i][1] >= staff_it_b[i-1][0]+in_trail:
        staff_it_b.append([eta_sort_b[i][1], eta_sort_b[i][2], eta_sort_b[i][3]])
    else:
        staff_it_b.append([staff_it_b[i-1][0]+in_trail, eta_sort_b[i][2], eta_sort_b[i][3]])
# print(eta_sort_b)
print('STAFF_it B:', staff_it_b)
rta_b = [[i[0]+transition_time, i[1], i[2], i[0]] for i in staff_it_b]
print('RTA B:', rta_b)

cp = []
for i in rta_a:
    cp.append(i)
for i in rta_b:
    cp.append(i)
cp = sorted(cp, key = lambda x: x[0])

cp = [[v[0], i, v[2], v[3]] for i, v in enumerate(cp)]

print(cp)
for i in range(0, len(staff_it_a)):
    plt.scatter(-1, eta_sort_a[i][1], marker = 's')
    plt.scatter(0, staff_it_a[i][0], marker = 's')

    plt.plot([0,1], [staff_it_a[i][0], rta_a[i][0]], linestyle = '--', color = 'k')
    plt.plot([-1, 0], [eta_sort_a[i][1], staff_it_a[i][0]], linestyle = '--', color = 'k')

    plt.text(-0.1, staff_it_a[i][0], f'{rta_a[i][2]}', weight = 'bold')
    plt.text(-1.1, eta_sort_a[i][1], f'{eta_sort_a[i][3]}', weight = 'bold')

for i in range(0, len(staff_it_b)):
    plt.scatter(2, staff_it_b[i][0], marker = 's')
    plt.scatter(3, eta_sort_b[i][1], marker = 's')

    plt.plot([2,1], [staff_it_b[i][0], rta_b[i][0]], linestyle = '--', color = 'k')
    plt.plot([2, 3], [staff_it_b[i][0], eta_sort_b[i][1]], linestyle = '--', color = 'k')

    plt.text(2, staff_it_b[i][0], f'{rta_b[i][2]}', weight = 'bold')
    plt.text(3.1, eta_sort_b[i][1], f'{eta_sort_b[i][3]}', weight = 'bold')

c = 0
for i in cp:
    print(f'{i[2]} ', end='')
    plt.scatter(1, i[0], marker = 's')
    if 'A' in i[2]:
        plt.text(0.9, i[0], f'{i[2]}', weight = 'bold')
    else:
        plt.text(1, i[0], f'{i[2]}', weight = 'bold')
    c+=1
plt.text(-1.1, 80, '$ETA_{FF}$\nGate A')
plt.text(-0.1, 80, '$STA_{FFit}$\nGate A')
plt.text(0.9, 80, '$RTA$\nRunway')
plt.text(1.9, 80, '$STA_{FFit}$\nGate B')
plt.text(2.9, 80, '$ETA_{FF}$\nGate B')
plt.xlim([-1.5, 3.5])
# for i in range(0, len(cp)):
#     if 'A' in cp[i][2]:
#         plt.plot([0, 1], [index_fcfs_b[i][4], index_eta_b[i][4]], color = 'black', linestyle = '--')
plt.show()


# %%
'''
CTAS Step 3 Computing STA
'''
sta = []
sta.append(cp[0])

for i in range(1, len(cp)):
    # print(cp[i][0], sta[i-1][0]+in_trail)
    if cp[i][0] > sta[i-1][0]+in_trail:
        print(cp[i][0], cp[i-1][0])
        sta.append(cp[i])
    else:
        sta.append([sta[i-1][0]+in_trail, cp[i][1], cp[i][2], cp[i][3]])
print(cp)
print(staff_it_a)
print(sta)

delay = []
print(sta[1])
for i in range(0, len(sta)):
    # print(sta[i][0], sta[i][3],sta[i][1], sta[i][2])
    delay.append([sta[i][0]-sta[i][3], sta[i][0], sta[i][1], sta[i][2], sta[i][3]])
print(delay)

# %%
'''
CTAS Step 4 Delay Distribution Function
'''
ddfC = []
ddfT = []
dtmax = 7
for i in delay:
    if i[0] <= dtmax:
        ddfC.append([0, i[1], i[2], i[3]])
        ddfT.append(i)
    else:
        ddfC.append([i[0]-dtmax, i[1], i[2], i[3], i[4]])
        ddfT.append([dtmax, i[1], i[2], i[3], i[4]])

# for i in range(0, len(delay)):
#     if ddfC[i][0]+ddfT[i][0] == delay[i][0]:
#         print('True!')
#     else:
#         print('False at index', i)
print(ddfC)
print(ddfT)

# %%
'''
CTAS Step 5 Scheduled Time of Arrival at Meter Gates
'''
staff_b = []
staff_a = []
# print(sta)
# print(ddfC)
sta_b = []
sta_a = []

for i in range(0, len(sta)):
    # print(sta[i][2], ddfC[i][3])
    if 'A' in sta[i][2]:
        sta_a.append([sta[i], ddfC[i][0]])
    else:
        sta_b.append([sta[i], ddfC[i][0]])
# print(sta_a)
# print(sta_b)

staff_b.append([sta_b[0][0][0]-transition_time, sta_b[0][0][1], sta_b[0][0][2]])
# print(staff_it_b[1], sta_b[1][1])
staff_b.append([staff_it_b[1][0]+sta_b[1][1], sta_b[1][0][1], sta_b[1][0][2]])
# print(staff_b)
for i in range(2, len(sta_b)):
    if staff_b[i-1][0]+in_trail > staff_it_b[i][0]+sta_b[i][1]:
        staff_b.append([staff_b[i-1][0]+in_trail, sta_b[i][0][1], sta_b[i][0][2]])
    else:
        staff_b.append([staff_it_b[i][0]+sta_b[i][1], sta_b[i][0][1], sta_b[i][0][2]])
# for i in sta_b:
print(staff_b)
staff_a.append([sta_a[0][0][0]-transition_time, sta_a[0][0][1], sta_a[0][0][2]])
# print(staff_it_a[1], sta_a[1][1])
staff_a.append([staff_it_a[1][0]+sta_a[1][1], sta_a[1][0][1], sta_a[1][0][2]])
# print(staff_a)
for i in range(2, len(sta_a)):
    if staff_a[i-1][0]+in_trail > staff_it_a[i][0]+sta_a[i][1]:
        staff_a.append([staff_a[i-1][0]+in_trail, sta_a[i][0][1], sta_a[i][0][2]])
    else:
        staff_a.append([staff_it_a[i][0]+sta_a[i][1], sta_a[i][0][1], sta_a[i][0][2]])
print(staff_a)

run = []
c = 0
j=0
while c == 0:
    # print(i,j
    try:
        if staff_b[c][0] < staff_a[c][0]:
            run.append(staff_b[c])
            staff_b.pop(0)
            j = 'b'
        else:
            run.append(staff_a[c])
            staff_a.pop(0)
            j = 'a'
    except:
        if j == 'b':
            run.append(staff_a[0])
        else:
            run.append(staff_b[0])
        c=1
    
print(run)
print(cp)
for i in cp:
    print(f'{i[2]} ', end='')
    plt.scatter(1, i[0], marker = 's')
    if 'A' in i[2]:
        plt.text(0.995, i[0], f'{i[2]}', weight = 'bold')
    else:
        plt.text(1, i[0], f'{i[2]}', weight = 'bold')
    c+=1
for i in run:
    print(f'{i[2]} ', end='')
    plt.scatter(1.05, i[0], marker = 's')
    if 'A' in i[2]:
        plt.text(1.045, i[0], f'{i[2]}', weight = 'bold')
    else:
        plt.text(1.05, i[0], f'{i[2]}', weight = 'bold')
    c+=1
for i in sta:
    print(f'{i[2]} ', end='')
    plt.scatter(1.025, i[0], marker = 's')
    if 'A' in i[2]:
        plt.text(1.02, i[0], f'{i[2]}', weight = 'bold')
    else:
        plt.text(1.025, i[0], f'{i[2]}', weight = 'bold')
    c+=1
plt.xlim(0.97, 1.08)
# plt.text(-1.1, 80, '$ETA_{FF}$\nGate A')
# plt.text(-0.1, 80, '$STA_{FFit}$\nGate A')
plt.text(0.995, 100, 'RTA\nRunway')
plt.text(1.02, 100, 'Calculated\nSTA')
plt.text(1.045, 100, 'Step 5\nRTA')
plt.show()

# %%



