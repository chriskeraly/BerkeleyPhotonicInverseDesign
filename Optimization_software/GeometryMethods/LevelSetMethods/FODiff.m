function [ Diff ] =FODiff(phi,direction,updown,dx)

Diff=zeros(size(phi));

if strcmp(updown,'down')+strcmp(updown,'up')+strcmp(updown,'central')==0
    a=[1,1]+[1,1,1];
end

if strcmp(direction,'y')
    if strcmp(updown,'up')
        Diff(1:end-1,:)=phi(2:end,:)-phi(1:(end-1),:);
        Diff(end,:)=0.0*(Diff(end-1,:));
        Diff=Diff/dx;
    else if strcmp(updown,'down')
            Diff(2:end,:)=phi(2:end,:)-phi(1:(end-1),:);
            Diff(1,:)=0.0*(Diff(2,:));
            Diff=Diff/dx;
        else
            Diff(2:end-1,:)=(phi(3:end,:)-phi(1:end-2,:))/2;
            Diff(1,:)=0.0*(phi(2,:)-phi(1,:));
            Diff(end,:)=0.0*(phi(end,:)-phi(end-1,:));
            Diff=Diff/dx;
        end
    end
    Diff(:,1)=0;
    Diff(:,end)=0;
else
    phi=phi';
    Diff=Diff';
    
    if strcmp(updown,'up')
        Diff(1:end-1,:)=phi(2:end,:)-phi(1:(end-1),:);
        Diff(end,:)=0.0*(Diff(end-1,:));
        Diff=Diff/dx;
    else if strcmp(updown,'down')
            Diff(2:end,:)=phi(2:end,:)-phi(1:(end-1),:);
            Diff(1,:)=0.0*(Diff(2,:));
            Diff=Diff/dx;
        else
            Diff(2:end-1,:)=(phi(3:end,:)-phi(1:end-2,:))/2;
            Diff(1,:)=phi(2,:)-phi(1,:);
            Diff(end,:)=0.0*(phi(end,:)-phi(end-1,:));
            Diff=Diff/dx;
        end
    end
    Diff(:,1)=0;
    Diff(:,end)=0;
    Diff=Diff';
end


end

