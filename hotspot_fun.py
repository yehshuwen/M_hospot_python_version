class hotspot_function(object):    
    
    def each_position(length, mut_data):
        tt = length+1
        TL = []
        TL = [x for x in range(1,tt)]
        #TL = np.zeros([tt,1],dtype=int)
        #def zxc(y):
        #num = len(where(mut_data['v2'] == y))
        #return(num)
        Count = np.zeros([tt,1],dtype=int)
        Count = pd.DataFrame(Count)
        Count.rename(columns={0: 'mutation'}, inplace=True)
        #i=1
        for i in range(len(mut_data)):
            if Count['mutation'][i]==0:
                num = np.size(np.where(mut_data['v2'] == mut_data['v2'].iloc[i]))
                #num = len(mut_data['v2'] == mut_data['v2'].iloc[i])
                position = mut_data['v2'].iloc[i]-1
                Count['mutation'][position] = num
                #Count.append(num)
                #Count = pd.DataFrame(Count)
                TL = pd.DataFrame(TL)
                rCount = pd.concat([TL,Count], axis=1, ignore_index=True)
                return rCount
    
    def Extend_region(peak_region_min, peak_region_max, remain_gap, mutation):
        if len(mutation)==1:
            peak_region = 10
        else:
            count = mutation.drop(mutation.index[-1])
            count1 = mutation.drop(mutation.index[0])
            #df.drop(df.index[2])
            gap1 = count1.values - count.values
            gap1 = list(chain.from_iterable(gap1))
            gap = sorted(gap1)
            new_gap = gap[0:math.ceil(len(gap)*remain_gap)]
            if len(new_gap) == 1:
                peak_region = math.ceil(mean(new_gap))
            else:
                peak_region = math.ceil(mean(new_gap)+std(new_gap))
        if peak_region < peak_region_min:
            peak_region = peak_region_min
        if peak_region > peak_region_max:
            peak_region = peak_region_max
        return peak_region
            
                            
    def Density_information(length, alph1, rCount):
        final = np.zeros([length,1],dtype=int)
        final = pd.DataFrame(final)
        #i = 1
        for i in range(length):
            final_score = np.zeros([1,1],dtype=int)
            final_score = pd.DataFrame(final_score)
            final_score.rename(columns={0: 'mutation'}, inplace=True)
            #k = 0
            for k in range(length):
                score = np.exp((-np.abs(i-k))/(alph1**2))*rCount.iloc[k:k+1,1:2]
                final_score = score.iloc[0:1,0:1] + final_score['mutation'].iloc[0]
            final.iloc[i:i+1,0:1]= final_score['mutation'].iloc[0]
        return final         
    
    
    def mountain_clusting(mutation, final, num_center, length, beta):
        beta_value = beta / (alph1 ** 2)
        final_p = []
        #final.index(final.values.max())
        P1 = final.score.idxmax()+1
        #w = 0
        PP = []
        P_score = []
        P_score = [x for x in range(length)]
        P_score = pd.DataFrame(P_score)
        P_score.rename(columns={0: 'score'}, inplace=True)
        for w in range(num_center):
            if w == 0:
                #PP[w] = P1
                PP.append(P1)
            else:
                #j = 0
                for j in range(length):
                    if(w == 1):
                        P_score['score'].iloc[j] = (final['score'].iloc[j] - ((final.values.max()) * (np.exp(-(beta_value * (P1 - j) ** 2)))))
                    else:
                        P_score['score'].iloc[j] = ( P_score['score'].iloc[j] - ((P_score.values.max()) * (np.exp(-(beta_value * (P - j) ** 2)))))
                P = P_score.score.idxmax()+1
                PP.append(P)
            final_p = np.unique(PP)
        
        final_p = mutation[mutation['position'].isin(final_p)]
       #df[df['A'].isin([3, 6])]
        return final_p

    def extend(final_P, peak_region, mutation, rCount):
        position_star = []
        position_end = []
        mutation1 = mutation
        mutation1.index = range(1,np.size(mutation1)+1)
        final_P.index = range(1,np.size(final_P)+1)
        e = 1
        for e in range(1,np.size(final_P)+1):
            grade = mutation['position']-final_P['position'].iloc[e-1]
            grade.index = range(1,np.size(grade)+1)
            center = grade[grade == 0].index[0]
            #center = np.where(grade == 0).index[0]
            if center == 1 or center == np.size(grade):
                if center == 1:
                    position_star.append(mutation['position'][center])
                    o=1
                    for o in range(1,np.size(mutation) - center+1):
                        if o == 1:
                            if np.abs(grade[center + o]) <= peak_region:
                                position_end.append(mutation['position'][center+o])
                            else:
                                position_end.append(mutation['position'][center])
                                break
                        else:
                            if np.abs(np.abs(grade[center+o]) - np.abs(grade[center + o -1])) <= peak_region:
                                position_end.append(mutation['position'][center+o])
                            else:
                                position_end.append(mutation['position'][center+o-1])
                                break
                else:
                    position_end.append(mutation['position'][center])
                    u = 1
                    for u in range(1,center-1+1):
                        if u == 1:
                            if np.abs(grade[center - u]) <= peak_region:
                                position_star.append(mutation['position'][center - u])
                            else:
                                position_star.append(mutation['position'][center])
                                break
                        else:
                            if np.abs(np.abs(grade[center - u]) - np.abs(grade[center - u +1])) <= peak_region:
                                position_star.append(mutation['position'][center - u])
                            else:
                                position_star.append(mutation['position'][center - u+1])
                                break
            else:
                u = 1
                for u in range(1,center-1+1):
                    if u == 1:
                        if np.abs(grade[center - u]) <= peak_region:
                            position_star.append(mutation['position'][center - u])
                        else:
                            position_star.append(mutation['position'][center])
                            break
                    else:
                        if np.abs(np.abs(grade[center - u]) - np.abs(grade[center -u +1])) <= peak_region:
                            position_star.append(mutation['position'][center - u])
                        else:
                            position_star.append(mutation['position'][center - u +1])
                            break
                    
                for o in range(1,np.size(mutation)-center+1):
                    if o == 1:
                        if np.abs(grade[center + o]) <= peak_region:
                            position_end.append(mutation['position'][center+o])
                        else:
                            position_end.append(mutation['position'][center])
                            break
                    else:
                        if np.abs(np.abs(grade[center + o]) - np.abs(grade[center + o -1])) <= peak_region:
                            position_end.append(mutation['position'][center+o])
                        else:     
                            position_end.append(mutation['position'][center+o-1])
                            break
        test_region = []
        g=0
        for g in range(np.size(position_star)):
            test_region.append(("%s_%s" % (int(position_star[g]), int(position_end[g]))))
        
        uniq_star = np.unique(position_star)
        uniq_end = np.unique(position_end)
        m_count = []
        for s in range(np.size(uniq_star)):
            #sum(rCount['mutation'].iloc[int(uniq_star[0])],rCount['mutation'].iloc[int(uniq_end[0])])
            if uniq_star[s] != uniq_end[s]:                                   
                mut = sum([rCount['mutation'][int(uniq_star[s])-1],rCount['mutation'][int(uniq_end[s])-1]])
                m_count.append(mut)
            else:
                mut = sum(rCount['mutation'][int(uniq_star[s])-1])
                m_count.append(mut)
        test_region = pd.DataFrame(test_region)
        m_count = pd.DataFrame(m_count)
        final_count = pd.concat([test_region,m_count], axis=1, ignore_index=True)                        
        fc4 = []
        for w in range(np.size(uniq_star)):
             leng = m_count[0].iloc[w]/(uniq_end[w] - uniq_star[w] + 1)
             fc4.append(leng)
        fc4 = pd.DataFrame(fc4)    
        final_count = pd.concat([final_count,fc4], axis=1, ignore_index=True)
        
        hotspot_rate = []
        for w in range(np.size(uniq_star)):
             rate = m_count[0].iloc[w]/sum([m_count[0]])
             hotspot_rate.append(rate)
        hotspot_rate = pd.DataFrame(hotspot_rate)    
        final_count = pd.concat([final_count,hotspot_rate], axis=1, ignore_index=True)
        return final_count