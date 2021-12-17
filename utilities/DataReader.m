
classdef DataReader < handle
methods
    function data = next(this)
        data = this.nextImpl;
    end
    function data = previous(this)
        data = this.previousImpl;
    end
    function data = current(this)
        data = this.currentImpl;
    end
end

methods(Abstract, Access=protected)
    data = nextImpl(this);
    data = previousImpl(this);
    data = currentImpl(this);
end
end



