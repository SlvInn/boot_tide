function bopt = cut_boot_blk_length(bopt)
% S. Innocenti. 2019/09 - modified 2021/06
% simulate the random block lengths or get a vector of (equal) fixed lengths

switch bopt.lBlk_dist
        
    case 'unif'
        assert(numel(bopt.lBlk)==2,'MBB Block length vector must be a 2 element vector (min-max block lengths for an Uniform distribution)')

        % min-max block length in hours
        minlBlk = bopt.lBlk(1)*24;
        maxlBlk = bopt.lBlk(2)*24;

        rng(bopt.seed(2))
        lBlkrand = randi([minlBlk maxlBlk],expNbl,n_boot);


    case 'pois' 

        assert(numel(bopt.lBlk)==1,'MBB Block length parameter (poisson distribution) must be a scalar (mean block lengths)')
        meanlBlk = bopt.lBlk;
        rng(bopt.seed(2))
        lBlkrand = 24.*poissrnd(meanlBlk-1, expNbl,n_boot)+1;

    case 'geom' 

        assert(numel(bopt.lBlk)==1,'MBB Block length parameter (geometric distribution) must be a scalar (1/mean block lengths)')
        assert(bopt.lBlk>0 && bopt.lBlk<1,'MBB Block length parameter must be a prob in (0,1)')

        meanlBlk = bopt.lBlk;
        rng(bopt.seed(2))
        lBlkrand = 24.*geornd(meanlBlk, expNbl,n_boot);

    case  'fix'

        assert(numel(bopt.lBlk)==1,'MBB Block length (fixed length) must be a scalar')
        assert(bopt.lBlk>=1 && bopt.lBlk<floor(ltin/24),'MBB Block length must be in [1,T]')
        
        meanlBlk = bopt.lBlk;
        lBlkrand = 24.*repmat(meanlBlk, expNbl,n_boot);
end