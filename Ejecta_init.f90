Module Ejecta_init

  use realkind
  use constants

  implicit none
  
  Contains

! In place of Module SLSN_data


    Subroutine Initialize_Ejecta()
      !!! Uses i_choose_ejecta from Ejecta_choice_param to choose the ejecta type; places it in the generic ejecta object Ejecta_obj
      use Ejecta_choice_param, ONLY: i_choose_ejecta
      use Ejecta_generic, ONLY: Ejecta_obj
      implicit none
      
        Ejecta_obj%i_choose_ejecta = i_choose_ejecta
           
        
        call Ejecta_obj%initialize()
        call Ejecta_obj%initialize_radiation()   !! added later (was in separate module) 
        
    End Subroutine Initialize_Ejecta
    

End Module Ejecta_init
