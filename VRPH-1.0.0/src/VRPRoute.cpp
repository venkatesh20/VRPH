////////////////////////////////////////////////////////////
//                                                        //
// This file is part of the VRPH software package for     //
// generating solutions to vehicle routing problems.      //
// VRPH was developed by Chris Groer (cgroer@gmail.com).  //
//                                                        //
// (c) Copyright 2010 Chris Groer.                        //
// All Rights Reserved.  VRPH is licensed under the       //
// Common Public License.  See LICENSE file for details.  //
//                                                        //
////////////////////////////////////////////////////////////

#include "VRPH.h"
#include "randvals.h"    // included for hashing


VRPRoute::VRPRoute()
{
    ///
    /// Default constructor for the VRPRoute.
    /// 

    this->name=NULL;
    this->ordering=NULL;
    this->x=NULL;
    this->y=NULL;

}

VRPRoute::VRPRoute(int n)
{
    ///
    /// Allocates memory for the VRPRoute fields large enough 
    /// for an n node problem.
    ///

    this->name=new char[8*n];// Used to add as column
    this->ordering=new int[n];
    this->x=new double[n];
    this->y=new double[n];

}


VRPRoute::~VRPRoute()
{
    ///
    /// Destructor for the VRPRoute.
    ///

    if(this->name)
        delete [] this->name;
    if(this->ordering)
        delete [] this->ordering;
    if(this->x)
        delete [] this->x;
    if(this->y)    
        delete [] this->y;

}


int VRPRoute::hash(int salt)
{
    /// 
    /// Computes a hash of the route, 
    /// returning an integer in the range [0,HASH_TABLE_SIZE-1].
    ///

    // Associate a random integer with each node in the route
    // XOR them together and take the low log_2(HASH_TABLE_SIZE) bits

    // Return 0 if # customers is 0 (the "NULL" route)
    if(this->num_customers==0)
        return 0;

    int i;
    int chunk=CHUNKSIZE;
    omp_lock_t lck1;
    omp_init_lock(&lck1);
    if(this->ordering[num_customers-1]<this->ordering[0])
    {
        fprintf(stderr,"Route has %d customers\n",this->num_customers);
        fprintf(stderr,"end<start! %d<%d\n",this->ordering[num_customers-1],this->ordering[0]);
        #pragma omp parallel for shared(chunk) schedule(dynamic, chunk) private(i)
        for(i=0;i<this->num_customers;i++) {
            omp_set_lock(&lck1);
            fprintf(stderr,"%d ",this->ordering[i]);
            omp_unset_lock(&lck1);
         }
        fprintf(stderr,"Length=%f;\nLoad=%d\nObj=%fStart=%d\nEnd=%d\n",this->length,this->load,this->obj_val,
            this->start,this->end);
    
        report_error("%s: Error in route hash\n",__FUNCTION__);
    }

    omp_destroy_lock(&lck1);
    int val = 0;
  
    #pragma omp parallel for shared(chunk) schedule(dynamic, chunk) private(i) reduction(^:val)
    for(i=0;i<this->num_customers; i++) 
        val ^= (randvals[(salt + VRPH_ABS(this->ordering[i])+
        VRPH_ABS(this->ordering[VRPH_MIN((this->num_customers)-1,i+1)]) )% NUM_RANDVALS]);

    #pragma omp parallel for shared(chunk) schedule(dynamic, chunk) private(i) reduction(+:val)
    for(i=0;i<this->num_customers; i++)
        val+=this->ordering[i];

    val=val&(HASH_TABLE_SIZE-1);
    // Get the right # of bits
    
    return val;

}


void VRPRoute::create_name()
{
    ///
    /// Creates a name for the route as a string. Format is
    /// hashval1_hashval2_ordering (delimited by _)
    ///

    // Make sure we have memory for the name
    if(!this->name)
    {
        fprintf(stderr,"No memory allocated for the route name!\n");
        exit(-1);
    }

    int i;
    char temp[100];

    int chunk=CHUNKSIZE;
    // Write the hash value    
    sprintf(this->name,"%d_%d_",this->hash(SALT_1),this->hash(SALT_2));

    // Write the ordering
    #pragma omp parallel for shared(chunk,temp) schedule(dynamic, chunk) private(i)
    for(i=0;i<this->num_customers-1;i++)
    {
        #pragma omp critical 
      {
        sprintf(temp,"%d_",this->ordering[i]);
        strcat(this->name,(const char *)temp);
      }
    }
    sprintf(temp,"%d",this->ordering[this->num_customers-1]);
    strcat(this->name,(const char *)temp);
   
    return;

}

VRPRouteWarehouse::VRPRouteWarehouse()
{
    ///
    /// Default constructor for the route warehouse.
    ///

    this->hash_table_size=0;
    this->num_unique_routes=0;

}

VRPRouteWarehouse::VRPRouteWarehouse(int h_size)
{
    ///
    /// Constructor for the rote warehouse with h_size entries in the hash
    /// table.  Best if h_size is a power of 2.
    ///

    // Allocate memory for the hash_table
    this->hash_table=new struct htable_entry[h_size];

    this->hash_table_size=h_size;
    this->num_unique_routes=0;  
    int chunk=CHUNKSIZE;  
    int i;
    #pragma omp parallel for shared(chunk) schedule(dynamic, chunk) private(i)
    for(i=0;i<h_size;i++)
        this->hash_table[i].num_vals=0;
        
}


VRPRouteWarehouse::~VRPRouteWarehouse()
{
    ///
    /// Destructor for the route warehouse.
    ///

    delete [] this->hash_table;

}

void VRPRouteWarehouse::liquidate()
{
    ///
    /// Clears the hash table, removing all routes from the WH.
    ///

    this->num_unique_routes=0;
    int chunk=CHUNKSIZE;
    int i;
    #pragma omp parallel for shared(chunk) schedule(dynamic, chunk) private(i)
    for(i=0;i<this->hash_table_size;i++)
        this->hash_table[i].num_vals=0;
}


void VRPRouteWarehouse::remove_route(int hash_val, int hash_val2)
{
    ///
    /// Removes a particular route from the hash table, using the
    /// second hash value to remove the correct entry if we have duplicates
    /// at the same location in the table.
    ///

    int i,j;

    if(this->hash_table[hash_val].num_vals==0)
    {
        fprintf(stderr,"Tried to remove empty hash_table entry!\n");
        fprintf(stderr,"hash_val=%d\n;hash_val2=%d\n",hash_val,hash_val2);
        report_error("%s: Error removing route from WH\n",__FUNCTION__);
    }

     int chunk=CHUNKSIZE;
    // Look for the second hash value
    bool flag=false;
    #pragma omp parallel for shared(chunk,hash_val2) schedule(dynamic, chunk) private(i,j)
    for(i=0;i<this->hash_table[hash_val].num_vals;i++)
    {
        if(( this->hash_table[hash_val].hash_val_2[i]==hash_val2) && (flag==false)) 
        {
            // This is the one to remove-shift everyone else down one
            // and decrement num_vals
            for(j=i+1;j<this->hash_table[hash_val].num_vals;j++)
            {
                // [j] -> [j-1]
                #pragma omp critical 
               { 
                this->hash_table[hash_val].length[j-1]=
                    this->hash_table[hash_val].length[j];
                this->hash_table[hash_val].hash_val_2[j-1]=
                    this->hash_table[hash_val].hash_val_2[j];
               }
            }

            this->hash_table[hash_val].num_vals--;

            this->num_unique_routes--;

            //return;
            flag=true;
            #pragma omp flush(flag)
        }
    }
    
    if (flag==true) return;

    // We shouldn't get here!
    // Annoying...
    fprintf(stderr,"Never found the route when removing from WH!!!\n");
    fprintf(stderr,"Looking for (%d, %d)\n",hash_val, hash_val2);
    for(i=0;i<this->hash_table[hash_val].num_vals;i++)
        fprintf(stderr,"Position %d:  (%d, %d, %f)\n",i, hash_val, 
            this->hash_table[hash_val].hash_val_2[i], hash_table[hash_val].length[i]);
    
    report_error("%s: Error removing route from WH\n",__FUNCTION__);

}


int VRPRouteWarehouse::add_route(VRPRoute *R)
{
    ///
    /// Adds the given route R to the warehouse, returning true
    /// if the addition is made and false otherwise.
    ///

    // Hash the route
    int i;
    R->hash_val  = R->hash(SALT_1);
    R->hash_val2 = R->hash(SALT_2);
    int hval=R->hash_val;
    int hval2=R->hash_val2;
    int chunk=CHUNKSIZE;
    omp_lock_t lck1;
    omp_init_lock(&lck1);
    
    bool flag1=false;
    bool flag2=false;

    if(this->hash_table[hval].num_vals>0)
    {
        // The entry is not vacant-must investigate further...
        // Make sure that we haven't filled up all the entries--should NEVER happen...
        // as long as the hash_table is big enough
        if(this->hash_table[hval].num_vals==NUM_ENTRIES)
        {
            fprintf(stderr, "Route hash table is too small!! \n");
            fprintf(stderr, "At entry %d\n",hval);
            fflush(stderr);
            #pragma omp parallel for shared(chunk,hval) schedule(dynamic, chunk) private(i)
            for(i=0;i<this->hash_table[hval].num_vals;i++)
            {
                omp_set_lock(&lck1);
                fprintf(stderr,"%d %f\n",i,this->hash_table[hval].length[i]);
                fflush(stderr);
                omp_unset_lock(&lck1);
            }
            fflush(stderr);
            report_error("%s: Error adding route to WH\n",__FUNCTION__);
        }        

        // Otherwise, there is room in the table
        for(i=0;i<this->hash_table[hval].num_vals;i++)
        {
            // Compare the second hash value-
            if ((this->hash_table[hval].hash_val_2[i] == hval2) && ((flag1==false)|| (flag2==false)))
            {
                
                // These surely must cover the same nodes - now see if we have 
                // a lower cost tour
                if(R->length<this->hash_table[hval].length[i] && 
                    VRPH_ABS(R->length - this->hash_table[hval].length[i])>VRPH_EPSILON)
                {
                    this->hash_table[hval].length[i]=R->length;
                    //return BETTER_ROUTE;
                    flag1=true;
                    #pragma omp flush(flag1)
                }
                else  {
                   // return DUPLICATE_ROUTE;
                   flag2=true;
                    #pragma omp flush(flag1)
                    // The column/route was the same or worse.
                }
                
            }
        }

      if ((flag1==true) && (flag2==false)) {
        return BETTER_ROUTE;
      } else if  ((flag1==false) && (flag2==true)) {
        return DUPLICATE_ROUTE;
      }
        // We didn't match any of the previous entries with this hval
        // New entry at this location
        this->hash_table[hval].length[this->hash_table[hval].num_vals]=
            R->length ;
        this->hash_table[hval].hash_val_2[this->hash_table[hval].num_vals]=
            hval2;
        
        this->hash_table[hval].num_vals++;
        this->num_unique_routes++;
        return ADDED_ROUTE;

    }
     omp_destroy_lock(&lck1);

    // ELSE...
    // This is a new entry--update the hash table
    this->hash_table[hval].num_vals=1;
    // Insert the length and h2
    this->hash_table[hval].length[0]=R->length ;
    this->hash_table[hval].hash_val_2[0]=hval2;
    this->num_unique_routes++;
    return ADDED_ROUTE;
}




