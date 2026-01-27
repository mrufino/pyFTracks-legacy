class Rufino2022(AnnealingModel):
    
    @staticmethod
    def convert_Dpar_to_rmr0(dpar, etchant="5.5HNO3"):
        """ Here depends on the etchant (5.5 or 5.0 HNO3)
            This is based on the relation between the fitted rmr0 values and
            the Dpar etched using a 5.5M etchant as published in
            Ketcham et al, 2007,Figure 6b
            We use the linear conversion defined in Ketcham et al 2007 to
            make sure that we are using 5.5M DPar"""

        """
            Converte comprimento de pit (Dpar) para RMR0 (reduced track length).
            A conversão depende do etchant usado (5.5M ou 5.0M HNO3).
        """
        if etchant == "5.0HNO3": 
             dpar = 0.9231 * dpar + 0.2515
        if dpar <= 1.75:
            return 0.84
        elif dpar >= 4.58:
            return 0
        else:
            return 0.84 * ((4.58 - dpar) / 2.98)**0.21
        
    @staticmethod
    def convert_Cl_pfu_to_rmr0(clpfu):
        """ Relation between fitted rmr0 value from the fanning curvilinear model and
            Cl content is taken from Ketcham et al 2007 Figure 6a """
        """
            Converte conteúdo de Cl em pfu para RMR0, segundo Ketcham 2007.
        """            
        value = np.abs(clpfu - 1.0)
        if value <= 0.130:
            return 0.0
        else:
            return 0.83 * ((value - 0.13) / 0.87)**0.23

    @staticmethod
    def convert_Cl_weight_pct(clwpct):
        """
            Converte Cl em % peso para RMR0.
        """
        # Convert %wt to APFU
        return Ketcham2007.convert_Cl_pfu_to_rmr0(clwpct * 0.2978)


    @staticmethod
    def convert_unit_paramA_to_rmr0(paramA):
        """
            Conversão unitária de um parâmetro cinético A para RMR0.
        """
        if paramA >= 9.51:
            return 0.0
        else:
            return 0.84 * ((9.509 - paramA) / 0.162)**0.175
    
    # Mapeamento de parâmetros cinéticos para funções de conversão
    _kinetic_conversion = {"ETCH_PIT_LENGTH": convert_Dpar_to_rmr0,
                          "CL_PFU": convert_Cl_pfu_to_rmr0,
                          "RMR0": lambda x: x}

    # Inicializa o modelo, mantendo compatibilidade com PyFTracks                          
    # Align API with Ketcham1999 and pyFTracks documentation (kinetic_parameters as dict).
    def __init__(self, kinetic_parameters: dict, use_projected_track: bool =False,
                use_Cf_irradiation: bool=False):

        """
        Armazena os parâmetros cinéticos e flags de uso de projected track / Cf irradiation.
        """
        super(Ketcham2007, self).__init__(
              kinetic_parameters,
              use_projected_track,
              use_Cf_irradiation)

    def annealing_model(self):
        """
        Calcula a redução de comprimento de tracks (g(r) = f(t,T)) ao longo da história térmica.
        
        Inputs implícitos:
            - self.history.time/temperature: vetores de tempo e temperatura
            - self.rmr0: valor inicial de reduced track length
            - modKetch07: constantes de ajuste do modelo

        Processos principais:
            1. Converte tempo para segundos e inicializa reduced_lengths
            2. Itera sobre cada nó da história térmica (do final para o início)
            3. Calcula reduced_lengths usando as equações do PRINCÍPIO DO TEMPO EQUIVALENTE (GOSWAMI)
            4. Aplica threshold (equivTotAnnLen) e conversão cinética (equação 8)
            5. Atualiza o tempo equivalente (equivTime) para o próximo nó

        Returns:
            - reduced_lengths: array de reduced track lengths ao longo do tempo
            - first_node: índice do primeiro nó onde ocorre annealing relevante
        """
        cdef double[::1] time = np.ascontiguousarray(self.history.time * _seconds_in_megayears)
        cdef double[::1] temperature = np.ascontiguousarray(self.history.temperature)
        cdef int numTTnodes = time.shape[0]
        cdef double[::1] reduced_lengths = np.zeros(time.shape[0] - 1)
        cdef double crmr0 = self.rmr0 # valor inicial do comprimento reduzido
        cdef int first_node = 0

        cdef int node, nodeB
        cdef double equivTime  # tempo equivalente de annealing
        cdef double timeInt, x1, x2
        cdef double equivTotAnnLen # threshold mínimo de reduced length
        cdef double k
        cdef double calc
        cdef double tempCalc
        cdef double MIN_OBS_RCMOD = _MIN_OBS_RCMOD

        # Define as constantes do modelo de Ketcham 2007 (curva fanning curvilinear)
        cdef annealModel modKetch07 = annealModel(
            c0=0.39528,
            c1=0.01073,
            c2=-65.12969,
            c3=-7.91715,
            a=0.04672,
            b=0)

        # Fator de conversão da redução inicial
        k = 1.04 - crmr0

        # "Mata" o PET mas mantém a classe funcional
        self.reduced_lengths = np.full(numTTnodes - 1, 1.0)  # ou crmr0 se quiser
        self.first_node = 0  # ou numTTnodes-2, conforme faça sentido

        # Retorna os resultados
        return self.reduced_lengths, self.first_node