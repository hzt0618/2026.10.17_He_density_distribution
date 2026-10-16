function lt=UtToLt(glong,ut)
lt=ut+glong/15;
lt(lt>24)=lt(lt>24)-24;
lt(lt<0)=lt(lt<0)+24;
end