classdef Adjust_Plots
    methods (Static)
        function Mean(f,obj)
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4 ];
            figure(f)
            set(f,'position',[237 1113 643 222])
            
            for j=1:2
                subplot(1,2,j)
                set(gca,'yscale','log','fontsize',14,'ylim',[1e-1,1e2]);
                pla = get(gca,'children');
                switch obj.Spatial
                    case {'NonSpatial'}
                        pla = get(gca,'children');
                        for i=1:3
                            set(pla(i),'marker','^','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(4),'color',COL(4,:),'linestyle','-','linewidth',4);
                        
                    case {'Spatial','Fixed'}
                        for i=1:3
                            set(pla((i-1)*2+1),'marker','^','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*2+2),'marker','o','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(7),'color',COL(4,:),'linestyle','-','linewidth',4);
                        set(pla(8),'color',COL(5,:),'linestyle','-','linewidth',4);
                end
            end
        end
        function f=Mean_Joint(f,obj)
            %%
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4 ];
            figure(f)
            
            for j=1:2
                subplot(1,2,j)
                TF(j).ax = gca;
            end
            
            f123 = figure(123);
            set(f123,'name','Means vs. Time');
            for j=1:2
                switch obj.Spatial
                    case {'NonSpatial'}
                        set(f,'position',[195 1026 560 207]);
                        subplot(1,2,j)
                        copyobj(allchild(TF(j).ax),gca)
                    case {'Spatial','Fixed'}
                        set(f,'position',[137 1113 643 222])
                        subplot(2,2,j)
                        copyobj(allchild(TF(j).ax),gca)
                        subplot(2,2,j+2)
                        copyobj(allchild(TF(j).ax),gca)
               end
            end
            
            switch obj.Spatial
                case {'NonSpatial'}
                    for j=1:2
                        subplot(1,2,j)
                        set(gca,'yscale','log','fontsize',14,'ylim',[1e-1,1e2]);
                        pla = get(gca,'children');
                        for i=1:3
                            set(pla((i-1)*2+1),'marker','^','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*2+2),'marker','o','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(7),'color',COL(4,:),'linestyle','-','linewidth',4);
                        set(pla(8),'color',COL(5,:),'linestyle','-','linewidth',4);
                        switch j
                            case 1
                                ylabel('average CTT1/STL1')
                                xlabel('time (min)')
                                title('0.2M NaCl')
                            case 2
                                xlabel('time (min)')
                                title('0.4M NaCl')
                        end
                        
                    end
                    
                case {'Spatial','Fixed'}
                    for j=1:4
                        subplot(2,2,j)
                        set(gca,'yscale','log','fontsize',14,'ylim',[1e-1,1e2]);
                        pla = get(gca,'children');
                        switch j
                            case {1,2} % STL1
                                set(pla([3,4,7,8,11,12]),'visible','off');
                                for i=1:2:6
                                    set(pla((i-1)*2+1),'marker','^','markersize',12,'color',COL(floor((i+1)/2),:),...
                                        'markerfacecolor',COL(floor((i+1)/2),:),'linestyle','none');
                                    set(pla((i-1)*2+2),'marker','o','markersize',12,'color',COL(floor((i+1)/2),:),...
                                        'markerfacecolor',COL(floor((i+1)/2),:),'linestyle','none');
                                end
                                set(pla(13),'color',COL(4,:),'linestyle','-','linewidth',4);
                                set(pla(14),'color',COL(5,:),'linestyle','-','linewidth',4);
                                set(pla([15,16]),'visible','off')
                            case {3,4} % CTT1
                                set(pla([1,2,5,6,9,10]),'visible','off');
                                for i=2:2:6
                                    set(pla((i-1)*2+1),'marker','^','markersize',12,'color',COL(floor((i+1)/2),:),...
                                        'markerfacecolor',COL(floor((i+1)/2),:),'linestyle','none');
                                    set(pla((i-1)*2+2),'marker','o','markersize',12,'color',COL(floor((i+1)/2),:),...
                                        'markerfacecolor',COL(floor((i+1)/2),:),'linestyle','none');
                                end
                                set(pla(15),'color',COL(4,:),'linestyle','-','linewidth',4);
                                set(pla(16),'color',COL(5,:),'linestyle','-','linewidth',4);
                                set(pla([13,14]),'visible','off')
                        end
                        switch j
                            case 1
                                ylabel('average CTT1 nuc/cyt')
                                title('0.2M NaCl')
                            case 2
                                title('0.4M NaCl')
                            case 3
                                ylabel('average STL1 nuc/cyt')
                                xlabel('time (min)')
                                title('0.2M NaCl')
                            case 4
                                xlabel('time (min)')
                                title('0.4M NaCl')
                        end
                        
                    end
                    
            end
            f=f123;
        end
        function Vars(f,obj)
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4;0.8 0.8 0.8 ];
            figure(f)
            set(f,'position',[238 814 641 224])
            
            for j=1:2
                subplot(1,2,j)
                set(gca,'yscale','log','fontsize',14,'ylim',[1e-2,1e8]);
                pla = get(gca,'children');
                switch obj.Spatial
                    case {'NonSpatial'}
                        pla = get(gca,'children');
                        for i=1:3
                            set(pla(i),'marker','^','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(4),'color',COL(4,:),'linestyle','-','linewidth',4);
                        
                    case {'Spatial','Fixed'}
                        for i=1:3
                            set(pla((i-1)*3+1),'marker','^','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*3+2),'marker','o','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*3+3),'marker','s','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(10),'color',COL(4,:),'linestyle','-','linewidth',4);
                        set(pla(11),'color',COL(5,:),'linestyle','-','linewidth',4);
                        set(pla(12),'color',COL(6,:),'linestyle','-','linewidth',4);
                end
            end
        end
        function f=Vars_Joint(f,obj)
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4;0.8 0.8 0.8 ];
            figure(f)
            set(f,'position',[238 814 641 224])
            
            for j=1:2
                subplot(1,2,j)
                TF(j).ax = gca;
            end
            
            f124 = figure(124);
            set(f124,'name','Variances vs. Time');
           
            switch obj.Spatial
                case 'NonSpatial'
                    obj_tmp = obj;
                    obj_tmp.Spatial = 'Spatial';
                    Adjust_Plots.Vars(f,obj);
                    
                case {'Spatial','Fixed'}
                    set(f124,'position',[440   162   559   636]);
                    for j=1:2
                        for i=0:2
                            subplot(3,2,j+2*i)
                            copyobj(allchild(TF(j).ax),gca)
                        end
                    end
                    
                    
                    for j=1:6
                        subplot(3,2,j)
                        set(gca,'yscale','log','fontsize',14,'ylim',[1e-4,1e4]);
                        pla = get(gca,'children');
                        
                        switch j
                            case {1,2} % STL1
                                set(pla([2,4,5,6,7,9,10]),'visible','off');
                                set(pla([2,4,5,6,7,9,10]+10),'visible','off');
                                set(pla([2,4,5,6,7,9,10]+20),'visible','off');
                                set(pla([2,4,5,6,7,9,10]+30),'visible','off');
                                for i=1:3
                                    set(pla(10*(i-1)+1),'marker','^','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+3),'marker','o','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+8),'marker','s','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                end
                                set(pla(31),'color',COL(4,:),'linestyle','-','linewidth',4);
                                set(pla(33),'color',COL(5,:),'linestyle','-','linewidth',4);
                                set(pla(38),'color',COL(6,:),'linestyle','-','linewidth',4);
                            case {3,4} % CTT1
                                set(pla([1 2 3 4 6 8 9]),'visible','off');
                                set(pla([1 2 3 4 6 8 9]+10),'visible','off');
                                set(pla([1 2 3 4 6 8 9]+20),'visible','off');
                                set(pla([1 2 3 4 6 8 9]+30),'visible','off');
                                for i=1:3
                                    set(pla(10*(i-1)+5),'marker','^','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+7),'marker','o','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+10),'marker','s','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                end
                                set(pla(35),'color',COL(4,:),'linestyle','-','linewidth',4);
                                set(pla(37),'color',COL(5,:),'linestyle','-','linewidth',4);
                                set(pla(40),'color',COL(6,:),'linestyle','-','linewidth',4);
                            case {5,6} % Both Genes
                                set(pla([1 3 5 7 8 10]),'visible','off');
                                set(pla([1 3 5 7 8 10]+10),'visible','off');
                                set(pla([1 3 5 7 8 10]+20),'visible','off');
                                set(pla([1 3 5 7 8 10]+30),'visible','off');
                                for i=1:3
                                    set(pla(10*(i-1)+2),'marker','^','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+4),'marker','o','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+6),'marker','s','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                    set(pla(10*(i-1)+9),'marker','x','markersize',8,'color',COL(i,:),...
                                        'markerfacecolor',COL(i,:),'linestyle','-','linewidth',1);
                                end
                                set(pla(32),'color',COL(6,:),'linestyle','-','linewidth',4);
                                set(pla(34),'color',COL(6,:),'linestyle','-','linewidth',4);
                                set(pla(36),'color',COL(6,:),'linestyle','-','linewidth',4);
                                set(pla(39),'color',COL(6,:),'linestyle','-','linewidth',4);
                                
                                
                        end
                        switch j
                            case 1
                                ylabel('(co)vars CTT1')
                                %             xlabel('time (min)')
                                title('0.2M NaCl')
                            case 2
                                %             ylabel('average CTT1')
                                %             xlabel('time (min)')
                                title('0.4M NaCl')
                            case 3
                                ylabel('(co)vars STL1')
                                xlabel('time (min)')
                                title('0.2M NaCl')
                            case 4
                                %             ylabel('average STL1')
                                xlabel('time (min)')
                                title('0.4M NaCl')
                            case 5
                                ylabel('co-vars (STL1/CTT1)')
                                xlabel('time (min)')
                                title('0.2M NaCl')
                            case 6
                                %             ylabel('average STL1')
                                xlabel('time (min)')
                                title('0.4M NaCl')
                        end
                        
                    end
                    f=f124;
            end
        end
        function On(f,spat)
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4; 0.8 0.8 0.8 ];
            figure(f)
            set(f,'position',[237 515 642 224])
            
            for j=1:2
                subplot(1,2,j)
                set(gca,'yscale','linear','fontsize',14,'ylim',[0,1.05]);
                pla = get(gca,'children');
                switch spat
                    case {'NonSpatial'}
                        set(gca,'yscale','linear','fontsize',14,'ylim',[0,1.05]);
                        pla = get(gca,'children');
                        for i=1:3
                            set(pla(i),'marker','^','markersize',12,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(4),'color',COL(4,:),'linestyle','-','linewidth',4);
                        
                    case {'Spatial','Fixed'}
                        set(gca,'yscale','log','fontsize',14,'ylim',[1e-4,1]);
                        for i=1:3
                            set(pla((i-1)*3+1),'marker','^','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*3+2),'marker','o','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                            set(pla((i-1)*3+3),'marker','s','markersize',8,'color',COL(i,:),...
                                'markerfacecolor',COL(i,:),'linestyle','none');
                        end
                        set(pla(10),'color',COL(4,:),'linestyle','-','linewidth',4);
                        set(pla(11),'color',COL(5,:),'linestyle','-','linewidth',4);
                        set(pla(12),'color',COL(6,:),'linestyle','-','linewidth',4);
                end
                switch j
                    case 1
                        ylabel('ON fraction')
                        xlabel('time (min)')
                    case 2
                        xlabel('time (min)')
                end

            end
        end
        function Marginal(f,obj,salt)
            COL = [0.4 0.4 0.8;0.8 0.4 0.8;0.4 0.8 0.8;0 0 0;0.4 0.4 0.4; 0.8 0.8 0.8 ];
            figure(f)
            set(f,'position',[175 915 560 420])
            
            for j=1:16
                subplot(4,4,j)
                set(gca,'xlim',[0 150],'ylim',[0 0.04],'yticklabel',[],'xticklabel',[],'fontsize',13);
                
                pla = get(gca,'children');
                for i=1:3
                    set(pla(i),'color',COL(i,:),'linewidth',2);
                end
                set(pla(4),'color',COL(4,:),'linestyle','-','linewidth',3);
                
            end
            %%
            for j=1:16
                subplot(4,4,j)
                text(75,0.035,obj.Gene);
                text(75,0.03,[num2str(salt),'M NaCl']);
                text(75,0.025,[num2str(obj.tt(j)),' min']);
                x = mod(j-1,4);y = floor((j-1)/4);
                set(gca,'Position',[0.01+x/4 0.77-y/4 0.23 0.23])
            end
        end
        function Joint(f,obj,salt,type,xax,yax)
            %%
            figure(f)
            STR = ['Joint ',xax,' vs ',yax,', ',num2str(salt),'M NaCl, ',type];
            set(f,'name',STR);
            set(f,'position',[173 392 560 420])
            %%
            for j=1:16
                subplot(4,4,j)
                set(gca,'yticklabel',[],'xticklabel',[],'fontsize',13);
                
                switch yax
                    case 'Nuc'
                        set(gca,'ylim',[0 10]);
                    case {'Cyt','CTT1','STL1'}
                        set(gca,'ylim',[0 200]);
                end
                switch xax
                    case 'Nuc'
                        set(gca,'xlim',[0 10]);
                    case {'Cyt','CTT1','STL1'}
                        set(gca,'xlim',[0 200]);
                end
                XLIM = get(gca,'xlim');
                YLIM = get(gca,'ylim');
                text(2/3*XLIM(2),7/8*YLIM(2),[num2str(obj.tt(j)),' min']);
                x = mod(j-1,4);y = floor((j-1)/4);
                set(gca,'Position',[0.01+x/4 0.77-y/4 0.23 0.23])
            end
        end
    end
end
